% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
trunkFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3.mat');
dendFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto-and-manual.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
asiRunId = '20190227T082543';

% For analysis of soma-distance
somaFile  = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_03_v2.mat');
wcFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto-and-manual.mat');

debugNmlDir = '';

% For the origin of this magic constant, see
% https://gitlab.mpcdf.mpg.de/connectomics/amotta/blob/534c026acd534957b84e395d697ac48b3cc6a7ad/matlab/+L4/+Spine/buildDendriteTrunkMask.m
trunkMinSize = 10 ^ 5.5;

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);

% NOTE(amotta): For translation of agglomerate → cell id
conn = load(connFile, 'axons', 'axonMeta', 'dendrites', 'denMeta', 'info');
conn.denMeta = connectEM.Dendrite.completeCellMeta(param, conn);

trunks = Util.load(trunkFile, 'dendrites');
trunks = Agglo.fromSuperAgglo(trunks);
trunkIds = Agglo.calculateVolume(param, trunks);
trunkIds = find(trunkIds > trunkMinSize);
trunks = trunks(trunkIds);

dendrites = load(dendFile);
shAgglos = dendrites.shAgglos;
dendrites = dendrites.dendrites;

% Only keep dendrites that also contain a trunk
curLUT = Agglo.buildLUT(maxSegId, trunks);
dendIds = Agglo.fromSuperAgglo(dendrites, true);
dendIds = find(cellfun(@(ids) any(curLUT(ids)), dendIds));
dendrites = dendrites(dendIds);

% Load ASI areas
[curDir, curAsiFile] = fileparts(connFile);
curAsiFile = sprintf('%s__%s_asiT.mat', curAsiFile, asiRunId);
curAsiFile = fullfile(curDir, curAsiFile);

asiT = Util.load(curAsiFile, 'asiT');
asiT = connectEM.Consistency.Calibration.apply(asiT);
asiT = asiT(asiT.area > 0, :);

% Sanity check
assert(numel(trunks) == numel(dendrites));

%% Map spine heads onto dendritic trunk
clear cur*;

shT = table;
shT.agglo = shAgglos;

[shT.length, shT.dendId, shT.trunkNodeId]  = ...
    connectEM.Dendrite.calculateSpineLengths( ...
        param, trunks, dendrites, shT.agglo);
    
% Sanity checks
assert(not(any(isnan(shT.length(shT.dendId > 0)))));
assert(not(any(isnan(shT.trunkNodeId(shT.dendId > 0)))));

% Map onto segment
curMask = not(isnan(shT.trunkNodeId));

shT.trunkSegId(:) = nan;
shT.trunkSegId(curMask) = arrayfun( ...
    @(d, n) dendrites(d).nodes(n, 4), ...
    shT.dendId(curMask), shT.trunkNodeId(curMask));

% Sanity check
assert(all(shT.trunkSegId(curMask)));

%% Calculate interspine distances
clear cur*;

curShT = shT;
curShT.id = reshape(1:height(curShT), [], 1);
curShT = curShT(curShT.dendId > 0, :);

curSynAgglos = repelem({[]}, height(shT), 1);
curSynAgglos(curShT.id) = num2cell(curShT.trunkSegId);

curAggloSynIds = accumarray( ...
    curShT.dendId, curShT.id, [numel(dendrites), 1], ...
    @(ids) {ids}, {zeros(0, 1, 'like', curShT.id)});

curAgglos = dendrites;

[dendShToShDists, dendShIds] = ...
    Synapse.calculateIntersynapseDistances( ...
        curAgglos, curAggloSynIds, curSynAgglos, ...
        'voxelSize', param.raw.voxelSize, ...
        'segWeights', ones(maxSegId, 1));

%% Debugging
clear cur*;
rng(0);

if ~isempty(debugNmlDir)
    curPoints = Seg.Global.getSegToPointMap(param);
    
    curDendIds = find(not(cellfun(@isempty, dendShToShDists)));
    curRandIds = curDendIds(randperm(numel(curDendIds)));

    curDendId = curRandIds(1);
    curDend = dendrites(curDendId);

    curShIds = dendShIds{curDendId};
    curShAgglos = shT.agglo(curShIds);
    curShTrunkNodeIds = shT.trunkNodeId(curShIds);

    curSkel = skeleton();
    curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    curSkel = Skeleton.setDescriptionFromRunInfo(curSkel, info);

    curComments = repelem({''}, size(curDend.nodes, 1), 1);
    curComments(curShTrunkNodeIds) = arrayfun( ...
        @(id) sprintf('Trunk to spine head %d', id), ...
        curShIds(:), 'UniformOutput', false);

    curSkel = curSkel.addTree( ...
        sprintf('Dendrite %d', curDendId), ...
        curDend.nodes(:, 1:3), curDend.edges, ...
        [], [], curComments);

    curShNodes = cellfun( ...
        @(ids) curPoints(ids, :), ...
        curShAgglos, 'UniformOutput', false);
    curShNames = arrayfun( ...
        @(id) sprintf('Spine head %d', id), ...
        curShIds, 'UniformOutput', false);

    curSkel = Skeleton.fromMST( ...
        curShNodes, param.raw.voxelSize, curSkel);
    curSkel.names(2:end) = curShNames;

    curSkel.write(fullfile(debugNmlDir, 'test.nml'));
end

%% Control: Compare ASI distributions across proximal dendrites
clear cur*;

curBw = 0.1;
curX = (-2):0.01:0.5;
curMinSynCount = 200;

curAsiT = asiT( ...
    asiT.type == 'PrimarySpine' ...
  & asiT.targetClass == 'ProximalDendrite', :);

curPdT = table;
[curPdT.aggloId, ~, curAsiT.postPdId] = unique(curAsiT.postAggloId);

curPdT.relAsiIds = accumarray( ...
    curAsiT.postPdId, ...
    reshape(1:height(curAsiT), [], 1), ...
    [], @(ids) {ids(:)}, {zeros(0, 1)});
curPdT.asiCount = cellfun( ...
    @numel, curPdT.relAsiIds);

% NOTE(amotta): Reduce sampling noise
curPdT = curPdT(curPdT.asiCount >= curMinSynCount, :);

% Plot results
curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');

% NOTE(amotta): Compute kernel density estimates of the probability density
% function, because visualization by histograms turned out to be difficult.
curData = nan(height(curPdT), numel(curX));
for curPdIdx = 1:height(curPdT)
    curPdAsiAreas = curPdT.relAsiIds{curPdIdx};
    curPdAsiAreas = curAsiT.area(curPdAsiAreas);
    curPdAsiAreas = log10(curPdAsiAreas);
    
    curY = ksdensity(curPdAsiAreas, curX, 'Bandwidth', curBw);
    curY = curY / sum(curY);
    
    curData(curPdIdx, :) = curY;
end

curColors = parula(height(curPdT));

% NOTE(amotta): Let's use the mode of the ASI area distribution to
% calculate the Z-index and the color of the proximal dendrite.
[~, curSortIds] = max(curData, [], 2);
[~, curSortIds] = sort(curX(curSortIds), 'ascend');

for curPdIdx = 1:height(curPdT)
    curPdId = curSortIds(curPdIdx);
    curColor = curColors(curPdIdx, :);
    
    plot(curAx, ...
        curX, curData(curPdId, :), ...
        'Color', curColor, 'LineWidth', 2);
end

xlabel(curAx, 'log10(ASI area [µm²])');
ylabel(curAx, 'Est. probability density');

curTitle = { ...
    'Primary spine synapses onto proximal dendrites'; sprintf( ...
    '(n = %d PDs with ≥ %d synapses)', height(curPdT), curMinSynCount)};
title(curAx, curTitle);

connectEM.Figure.config(curFig, info);
curFig.Position(3:4) = [390, 275];


% NOTE(amotta): Some more information on the most extreme neurons
curPdT.cellId = conn.denMeta.cellId(curPdT.aggloId);

fprintf([ ...
    '\n', ...
    'ASI area distribution over proximal dendrites\n', ...
    '* Neurons with smallest ASI area modes: %d, %d, %d\n', ...
    '* Neurons with largest ASI area modes: %d, %d, %d\n'], ...
    curPdT.cellId(curSortIds(1:3)), ...
    curPdT.cellId(curSortIds((end - 2):end)));

%% Let's have a look at synapses onto proximal dendrites
clear cur*;
rng(0);

% NOTE(amotta): If set, the analysis will be performed exclusively for the
% synapses onto this particular cell.
%   Most importantly, this also affects the null model: The synapse size
% distributions in the surround will be compared against other synapses
% onto the same neuron.
cellId = [];

% NOTE(amotta): The first run corresponds to actually observed data. All
% subsequent runs use randomly redistributed ASI areas.
numRuns = 1 + 10;
distThresh = 2500;

curShT = shT;
curAsiT = asiT( ...
    asiT.type == 'PrimarySpine' ...
  & asiT.targetClass == 'ProximalDendrite' ...
  & ismember(asiT.axonClass, {'Corticocortical', 'Thalamocortical'}), :);

% IMPORTANT(amotta): Use `postAggloId` before transforming it!
curAsiT.postCellId = conn.denMeta.cellId(curAsiT.postAggloId);
curAsiT.postAggloId = shT.dendId(curAsiT.shId);

if ~isempty(cellId)
    curAsiT = curAsiT(curAsiT.postCellId == cellId, :);
end

seedAreas = nan([height(curAsiT), numRuns]);
condAreas = cell([height(curAsiT), numRuns]);

for curRunIdx = 1:numRuns
    curShT.asiArea(:) = nan;
    curShT.asiArea(curAsiT.shId) = curAsiT.area;
    
    for curIdx = 1:height(curAsiT)
        curSeedShId = curAsiT.shId(curIdx);
        curSeedDendId = curAsiT.postAggloId(curIdx);
        
        curSeedAsiArea = curShT.asiArea(curSeedShId);

        % Find other spine heads onto same dendrite
        curDendShIds = dendShIds{curSeedDendId};
        curSeedShMask = curDendShIds == curSeedShId;
        assert(sum(curSeedShMask) == 1);

        % Restrict to other spine heads in surround
        curDendDists = dendShToShDists{curSeedDendId};
        curDendDists = curDendDists(:, curSeedShMask);

        curOtherShMask = curDendDists < distThresh;
        curOtherShMask = curOtherShMask & not(curSeedShMask);

        curOtherShIds = curDendShIds(curOtherShMask);
        curOtherAsiAreas = curShT.asiArea(curOtherShIds);
        curOtherAsiAreas = rmmissing(curOtherAsiAreas);
        
        seedAreas(curIdx, curRunIdx) = curSeedAsiArea;
        condAreas{curIdx, curRunIdx} = curOtherAsiAreas(:);
    end
    
    % Shuffle areas
    curRandIds = randperm(height(curAsiT));
    curAsiT.area = curAsiT.area(curRandIds);
end

%% Plot size versus number of neighboring spines
clear cur*;

curX = log10(seedAreas(:, 1));
curY = condAreas(:, 1);
curY = cellfun(@numel, curY);

curFit = fitlm(curX, curY);

curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');

scatter(curAx, curX, curY, '.');
curFitX = curAx.XLim(:);
curFitY = curFit.predict(curFitX);
plot(curAx, curFitX, curFitY, 'k--');
connectEM.Figure.config(curFig, info);

%% Plot results
clear cur*;

% Configuration
curImSize = 256;
curLims = [-1.5, 0];

curTicks = [-1, 0];
curTickLabels = {'-1', '0'};
curTickIds = linspace(curLims(1), curLims(2), curImSize);
[~, curTickIds] = arrayfun(@(t) min(abs(t - curTickIds)), curTicks);

curXLabel = 'log10(reference ASI area [µm²])';
curYLabel = 'log10(ASI area within %g µm [µm³])';
curYLabel = sprintf(curYLabel, distThresh / 1E3);

curRedCmap = connectEM.Figure.whiteRed(256);
curBlueRedCmap = connectEM.Figure.blueWhiteRed(256);

curTitle = '';
if ~isempty(cellId)
    curTitle = sprintf('Synapses onto neuron %d', cellId);
end

curConfigAxis = @(ax) set(ax, ...
    'PlotBoxAspectRatio', [1, 1, 1], ...
    'DataAspectRatioMode', 'auto', ...
    'XLim', [0, curImSize] + 0.5, 'XDir', 'normal', ...
    'XTick', curTickIds, 'XTickLabel', curTickLabels, ...
    'YLim', [0, curImSize] + 0.5, 'YDir', 'normal', ...
    'YTick', curTickIds, 'YTickLabel', curTickLabels);

curBw = [];
curDens = cell(numRuns, 1);
for curRunIdx = 1:numRuns
    curSeedAreas = seedAreas(:, curRunIdx);
    curCondAreas = condAreas(:, curRunIdx);
    
    curN = cellfun(@numel, curCondAreas);
    curW = repelem(1 ./ curN, curN);
    curX = repelem(curSeedAreas, curN);
    curY = cell2mat(curCondAreas);

    curX = log10(curX);
    curY = log10(curY);
    
   [curDens{curRunIdx}, curRunBw] = kde2d( ...
        curX, curY, curW, curImSize, curLims, curLims, curBw);
    
    % NOTE(amotta): The first run contains the "observed data". We use an
    % empty bandwidth vector for this run in order for `kde2d` to derive
    % these values.
    %   But for all subsequent runs we want to use the bandwidth derived
    % from the actualy data. So, we set it here.
    if isempty(curBw); curBw = curRunBw; end
    
    curDens{curRunIdx} = ...
        curDens{curRunIdx} ...
      / sum(curDens{curRunIdx}(:));
end

curRealDens = curDens{1};
curCtrlDens = mean(cat(3, curDens{2:end}), 3);
curDiffDens = curRealDens - curCtrlDens;

curRealCond = sum(curRealDens, 1);
curRealCond(~curRealCond) = eps;
curRealCond = curRealDens ./ curRealCond;

curCtrlCond = sum(curCtrlDens, 1);
curCtrlCond(~curCtrlCond) = eps;
curCtrlCond = curCtrlDens ./ curCtrlCond;
curDiffCond = curRealCond - curCtrlCond;


% Plot 1
curFig = figure();
curAx = axes(curFig);
imagesc(curAx, curRealDens);
curCbar = colorbar(curAx);

xlabel(curAx, curXLabel);
ylabel(curAx, curYLabel);
curAx.CLim = [0, max(curRealDens(:))];
curCbar.Label.String = 'P(ref. ASI, close-by ASIs)';

curConfigAxis(curAx);
colormap(curAx, curRedCmap);
if ~isempty(curTitle); title(curAx, curTitle); end
connectEM.Figure.config(curFig, info);


% Plot 2
curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');

imagesc(curAx, curDiffDens);
curCbar = colorbar(curAx);

xlabel(curAx, curXLabel);
ylabel(curAx, curYLabel);
curAx.CLim = [-1, +1] * max(abs(curDiffDens(:)));
curCbar.Label.String = 'ΔP(ref. ASI, close-by ASIs) to random';

curConfigAxis(curAx);
colormap(curAx, curBlueRedCmap);
plot(curAx, xlim(curAx), ylim(curAx), 'k--');
if ~isempty(curTitle); title(curAx, curTitle); end
connectEM.Figure.config(curFig, info);


% Plot 3
curFig = figure();
curAx = axes(curFig);
imagesc(curAx, curRealCond);
curCbar = colorbar(curAx);

xlabel(curAx, curXLabel);
ylabel(curAx, curYLabel);
curAx.CLim = [0, max(curRealCond(:))];
curCbar.Label.String = 'P(close-by ASIs | ref. ASI)';

curConfigAxis(curAx);
colormap(curAx, curRedCmap);
if ~isempty(curTitle); title(curAx, curTitle); end
connectEM.Figure.config(curFig, info);


% Plot 4
curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');
imagesc(curAx, curDiffCond);
curCbar = colorbar(curAx);

curClim = curDiffCond(:, 75:225);
curClim = [-1, +1] * max(abs(curClim(:)));

xlabel(curAx, curXLabel);
ylabel(curAx, curYLabel);
curAx.CLim = curClim;
curCbar.Label.String = 'ΔP(close-by ASIs | ref. ASI) to random';

curConfigAxis(curAx);
colormap(curAx, curBlueRedCmap);
plot(curAx, xlim(curAx), ylim(curAx), 'k--');
if ~isempty(curTitle); title(curAx, curTitle); end
connectEM.Figure.config(curFig, info);

%% Find size range with signs of homeostatic suppression in surround
clear cur*;

curSeedAreas = seedAreas(:, 1);
curCondAreas = condAreas(:, 1);

% To log10-space
curSeedAreas = log10(curSeedAreas);
curCondAreas = cellfun(@log10, curCondAreas, 'UniformOutput', false);

curCondN = cellfun(@numel, curCondAreas);
curCondAreas = cell2mat(curCondAreas);

curYs = [];
for curX = -1.5:0.05:0
    % NOTE(amotta): Uniform kernel
    curCondW = curSeedAreas - curX;
    curCondW = double(abs(curCondW) < 0.1);

    assert(isequal(size(curCondW), size(curCondN)));
    curCondW = repelem(curCondW ./ curCondN, curCondN);
    curCondW = curCondW / sum(curCondW(:));

    curY = ...
        sum(curCondAreas(:) .* curCondW(:)) ...
      - mean(curSeedAreas);
    curYs(end + 1) = curY;
end

figure;
hold on;
plot(curYs);
plot(xlim(), [0, 0], 'k--');

%% Look at effect of soma distance
clear cur*;

curWcs = load(wcFile, 'dendrites', 'indWholeCells');
curWcs = curWcs.dendrites(curWcs.indWholeCells);

curSomata = load(somaFile, 'dendAgglos', 'indSomata');
curSomata = curSomata.dendAgglos(curSomata.indSomata);

% Build soma → whole cell mapping
curWcLUT = Agglo.fromSuperAgglo(curWcs);
curWcLUT = Agglo.buildLUT(maxSegId, curWcLUT);
curSomaToWc = cellfun(@(ids) mode(nonzeros(curWcLUT(ids))), curSomata);
assert(isequaln(sort(curSomaToWc), unique(curSomaToWc)));

% Build whole cell → soma mapping
[~, curWcToSoma] = ismember(1:numel(curWcs), curSomaToWc);
assert(all(curWcToSoma));

% Calculate soma distances
somaDists = nan(maxSegId, 1);

for curIdx = 1:numel(curWcs)
    curWc = curWcs(curIdx);
    curSomaId = curWcToSoma(curIdx);
    curSoma = curSomata{curSomaId};
    
    curSomaNodeIds = ismember(curWc.nodes(:, 4), curSoma);
    curSomaNodeIds = find(curSomaNodeIds);
    
    curGraph = ...
        curWc.nodes(curWc.edges(:, 1), 1:3) ...
      - curWc.nodes(curWc.edges(:, 2), 1:3);
    curGraph = curGraph .* param.raw.voxelSize;
    curGraph = sqrt(sum(curGraph .* curGraph, 2));
    
    curGraph = graph( ...
        curWc.edges(:, 1), curWc.edges(:, 2), ...
        curGraph, size(curWc.nodes, 1));
    
    curSomaDists = distances(curGraph, curSomaNodeIds);
    curSomaDists = min(curSomaDists, [], 1);
    
    curSegIds = curWc.nodes(:, 4);
    curMask = not(isnan(curSegIds));
    
    curSegIds = curSegIds(curMask);
    curSomaDists = curSomaDists(curMask);
    
    somaDists(curSegIds) = curSomaDists;
    curIdx %#ok
end

%% Plot results
curAsiT = asiT( ...
    asiT.type == 'PrimarySpine' ...
  & asiT.targetClass == 'ProximalDendrite' ...
  & ismember(asiT.axonClass, {'Corticocortical', 'Thalamocortical'}), :);

curShT = shT;
curShT.asiArea(:) = nan;
curShT.asiArea(curAsiT.shId) = curAsiT.area;

curShT.somaDist(:) = nan;
curMask = ~isnan(curShT.trunkSegId);
curShT.somaDist(curMask) = somaDists(curShT.trunkSegId(curMask));

curShT = curShT( ...
    not(isnan(curShT.asiArea)) ...
  & not(isnan(curShT.somaDist)), :);

curX = curShT.somaDist / 1E3;
curY = log10(curShT.asiArea);

curFit = fitlm(curX(:), curY(:));
disp(curFit);


% Plot 1
curFig = figure;
curAx = axes(curFig);
hold(curAx, 'on');

scatter(curAx, curX, curY, '.');
curFitX = curAx.XLim(:);
curFitY = curFit.predict(curFitX);
plot(curAx, curFitX, curFitY, 'k--');

axis(curAx, 'square');
xlabel(curAx, 'Distance to soma [µm]');
ylabel(curAx, 'log10(ASI area [µm²]');
connectEM.Figure.config(curFig, info);

%% Utilities
function [im, bwOut] = kde2d(x, y, w, imSize, xlim, ylim, bwIn)
   [gridY, gridX] = ndgrid( ...
        linspace(ylim(1), ylim(2), imSize), ...
        linspace(xlim(1), xlim(2), imSize));
    grid = horzcat(gridY(:), gridX(:));

    optKvPairs = cell(0, 2);
    if ~isempty(bwIn); optKvPairs(end + 1, :) = {'Bandwidth', bwIn}; end
    if ~isempty(w); optKvPairs(end + 1, :) = {'Weights', w}; end
    optKvPairs = transpose(optKvPairs);
    
    im = cat(2, y(:), x(:));
   [im, ~, bw] = ksdensity(im, grid, optKvPairs{:});
    im = reshape(im, [imSize, imSize]);
        
    bwOut = bwIn;
    if isempty(bwOut); bwOut = bw; end
end
