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
curTrunkMask = Agglo.calculateVolume(param, trunks);
curTrunkMask = curTrunkMask > trunkMinSize;
trunks = trunks(curTrunkMask);

shT = table;
% NOTE(amotta): The dendrites of the connectome only consist of segment
% equivalence classes. But in this script, we're interested in the
% intersynapse and interspine distances along the dendritic trunk.
%   That's why we use the dendritic super-agglomerates and the spine head
% agglomerates from a separate file here.
[dendrites, shT.agglo] = Util.load( ...
    dendFile, 'dendrites', 'shAgglos');

% Only keep dendrites that also contain a trunk
curTrunkLUT = Agglo.buildLUT(maxSegId, trunks);
curDendMask = Agglo.fromSuperAgglo(dendrites, true);
curDendMask = cellfun(@(ids) any(curTrunkLUT(ids)), curDendMask);
dendrites = dendrites(curDendMask);

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

%% Clean up spine head table
% NOTE(amotta): As mentioned in the first sections of this script, we
% required dendrite super-agglomerates from a separate file to compute the
% intersynapse and interspine distances.
%   But now that we've done that, let's translate all dendrite indices and
% get rid of data that refers to dendrites other than the ones in the
% connectome file.
clear cur*;

curConnLUT = Agglo.buildLUT(maxSegId, conn.dendrites);
curConnIds = Agglo.fromSuperAgglo(dendrites);

curConnIds = cellfun(@(ids) mode(nonzeros(curConnLUT(ids))), curConnIds);
% NOTE(amotta): MATLAB's `mode` function returns NaN on empty inputs.
curConnIds(isnan(curConnIds)) = 0;

% Update dendrite IDs
curMask = shT.dendId > 0;
shT.dendId(curMask) = curConnIds(shT.dendId(curMask));

dendT = table;
dendT.shIds = repelem({zeros(0, 1)}, numel(conn.dendrites), 1);
dendT.shIds(curConnIds(curConnIds > 0)) = dendShIds(curConnIds > 0);
dendT.shToShDists = repelem({zeros(0, 0)}, numel(conn.dendrites), 1);
dendT.shToShDists(curConnIds(curConnIds > 0)) = dendShToShDists(curConnIds > 0);

% NOTE(amotta): Remove obsolete data
shT(:, {'trunkNodeId', 'trunkSegId'}) = [];
clear dendrites trunks dendShIds dendShToShDists;

%% Complete spine head table
clear cur*;

curSegSizes = Seg.Global.getSegToSizeMap(param);
shT.vol = cellfun(@(ids) sum(curSegSizes(ids)), shT.agglo);
shT.vol = shT.vol * prod(param.raw.voxelSize / 1E3);

curMask = shT.dendId ~= 0;
shT.targetClass(:) = categorical({'Unknown'});
shT.targetClass(curMask) = conn.denMeta.targetClass(shT.dendId(curMask));

%% Control: Compare ASI distributions across proximal dendrites
clear cur*;

bw = 0.1;
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
    
    curY = ksdensity(curPdAsiAreas, curX, 'Bandwidth', bw);
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

% Configurations
sizeVar = 'vol';
imSize = 256;

switch lower(sizeVar)
    case 'asi'
        lims = [-1.5, 0];
        bw = [0.1, 0.1];
        sizeLabel = 'ASI area [µm²]';
    case 'vol'
        lims = [-2.5, -0];
        bw = [0.2, 0.2];
        sizeLabel = 'spine head volume [µm³]';
    otherwise
        error('Unknown variable');
end

% NOTE(amotta): The first run corresponds to actually observed data. All
% subsequent runs use randomly redistributed ASI areas.
numRuns = 2;
distThresh = 2500;
minShCount = 0;

curShT = shT;
curShT.size = curShT.(sizeVar);
curShT.id = reshape(1:height(curShT), [], 1);

dendIds = curShT.targetClass == 'WholeCell';
dendIds = curShT.dendId(dendIds);

[dendIds, ~, curCount] = unique(dendIds);
curCount = accumarray(curCount, 1);
dendIds = dendIds(curCount > minShCount);

dendIds = setdiff(dendIds, 0);

tic;
dendData = cell(numel(dendIds), 1);
for curDendIdx = 1:numel(dendIds)
    curDendId = dendIds(curDendIdx);
    curDendShIds = dendT.shIds{curDendId};
    
    curDendMask = curShT.dendId == curDendId;
    curDendShT = curShT(curDendMask, :);
    
    % NOTE(amotta): In control runs, the spine sizes will be shuffled. But
    % that shuffling only affects `curShT`. To make sure that the spine
    % size will always be pulled from from there and not `curDendShT`, it
    % is probably safer to just remove it here.
    curDendShT.size = [];

    curDendSeedSizes = nan([height(curDendShT), numRuns]);
    curDendCondSizes = cell([height(curDendShT), numRuns]);

    rng(0);
    for curRunIdx = 1:numRuns
        for curIdx = 1:height(curDendShT)
            curSeedShId = curDendShT.id(curIdx);
            curSeedSize = curShT.size(curSeedShId);

            % Find other spine heads onto same dendrite
            curSeedShMask = curDendShIds == curSeedShId;
            assert(sum(curSeedShMask) == 1);

            % Restrict to other spine heads in surround
            curDendDists = dendT.shToShDists{curDendId};
            curDendDists = curDendDists(:, curSeedShMask);

            curOtherShMask = curDendDists < distThresh;
            curOtherShMask = curOtherShMask & not(curSeedShMask);

            curOtherShIds = curDendShIds(curOtherShMask);
            curOtherSizes = curShT.size(curOtherShIds);
            curOtherSizes = rmmissing(curOtherSizes);

            curDendSeedSizes(curIdx, curRunIdx) = curSeedSize;
            curDendCondSizes{curIdx, curRunIdx} = curOtherSizes(:);
        end

        % Shuffle sizes
        curShuffled = curShT.size(curDendMask);
        curShuffled = curShuffled(randperm(numel(curShuffled)));
        curShT.size(curDendMask) = curShuffled;
        clear curShuffled;
    end

    curDens = cell(numRuns, 1);
    for curRunIdx = 1:numRuns
        curSeedSizes = curDendSeedSizes(:, curRunIdx);
        curCondSizes = curDendCondSizes(:, curRunIdx);

        curN = cellfun(@numel, curCondSizes);
        curX = repelem(curSeedSizes, curN);
        curY = cell2mat(curCondSizes);

        % NOTE(amotta): As of now (27.01.2020), it's not yet clear to
        % me which of the following strategies is what we want to do:
        %
        % * Assign each synapse in each neighborhood the same weight?
        %   This would mean that one surround with many small synapses
        %   would dominate over another surround with one large synapse.
        % * Or, assign each synapse in each neighborhood that is
        %   inversely proportional to the number of synapses in the
        %   neighborhood? This would allow us to ask: What is the
        %   average size distribution in the surround of a synapse of
        %   size X?
        %
        % The latter seems closer to what I'm interested in. This
        % approach is more robust against the synapse density being
        % correlated with the size of the surround seed synapse.
        curW = [];

        curX = log10(curX);
        curY = log10(curY);

       [curDens{curRunIdx}, curRunBw] = kde2d( ...
            curX, curY, curW, imSize, lims, lims, bw);

        curDens{curRunIdx} = ...
            curDens{curRunIdx} ...
          / sum(curDens{curRunIdx}(:));
    end

    curDensReal = curDens{1};
    curDensCtrl = cat(3, curDens{2:end});
    curDensCtrl = mean(curDensCtrl, 3);

    curCellData = struct;
    curCellData.seedSizes = curDendSeedSizes(:, 1);
    curCellData.condSizes = curDendCondSizes(:, 1);
    curCellData.densReal = curDensReal;
    curCellData.densCtrl = curDensCtrl;

    dendData{curDendIdx} = curCellData;
    Util.progressBar(curDendIdx, numel(dendIds));
end

dendData = cat(1, dendData{:});


%% Plot results
clear cur*;

% NOTE(amotta): Set to `ctrlDens` for control runs
curPlotVar = 'densReal';

curTicks = [-2, -1, 0];
curTickLabels = {'-2', '-1', '0'};
curTickIds = linspace(lims(1), lims(2), imSize);
[~, curTickIds] = arrayfun(@(t) min(abs(t - curTickIds)), curTicks);

curXLabel = sprintf('log10(reference %s)', sizeLabel);
curYLabel = sprintf('log10(%s within %g µm)', sizeLabel, distThresh / 1E3);
curBlueRedCmap = connectEM.Figure.blueWhiteRed(256);

curConfigAxis = @(ax) set(ax, ...
    'PlotBoxAspectRatio', [1, 1, 1], ...
    'DataAspectRatioMode', 'auto', ...
    'XLim', [0, imSize] + 0.5, 'XDir', 'normal', ...
    'XTick', curTickIds, 'XTickLabel', curTickLabels, ...
    'YLim', [0, imSize] + 0.5, 'YDir', 'normal', ...
    'YTick', curTickIds, 'YTickLabel', curTickLabels);


curSpineCounts = nan(numel(dendData), 1);
curSpineDensDiffs = nan(numel(dendData), imSize);

curWeights = nan(numel(dendData), 1);
curDensDiffs = nan(numel(dendData), imSize, imSize);
curCondDiffs = nan(numel(dendData), imSize, imSize);
curCondDiffDiags = nan(numel(dendData), imSize);

for curDendIdx = 1:numel(dendData)
    curDendId = dendIds(curDendIdx);
    curCellData = dendData(curDendIdx);
    
    curWeight = sum(not(cellfun( ...
        @isempty, curCellData.condSizes)));
    curWeights(curDendIdx) = curWeight;
    if ~curWeight; continue; end
    
    curDensReal = curCellData.(curPlotVar);
    assert(not(any(isnan(curDensReal(:)))));
    assert(abs(sum(curDensReal(:)) - 1) < 1E-6);
    
    
    % Size distribution
    curDensSeed = ksdensity( ...
        log10(curCellData.seedSizes), ...
        linspace(lims(1), lims(2), imSize), ...
        'Bandwidth', bw(1));
    
    curDensSeed = curDensSeed / sum(curDensSeed(:));
    curDensSeed = reshape(curDensSeed, 1, []);
    
    
    % Check size-dependence of spine density
    curSpineCount = numel(curCellData.seedSizes);
    curSpineCounts(curDendIdx) = curSpineCount;
    
    curTotalSpinesInSurround = sum(cellfun(@numel, curCellData.condSizes));
    curAvgSpinesInSurround = curTotalSpinesInSurround / curSpineCount;
    
    curSpineDensCtrl = curAvgSpinesInSurround * curSpineCount;
    curSpineDensCtrl = curDensSeed(:) .* curDensSeed .* curSpineDensCtrl;
    assert(abs(sum(curSpineDensCtrl(:)) - curTotalSpinesInSurround) < 1);
    
    curSpineDensReal = curTotalSpinesInSurround * curDensReal;
    curSpineDensDiff = curSpineDensReal - curSpineDensCtrl;
    
    curSpineDensDiffs(curDendIdx, :) = ...
        sum(curSpineDensDiff, 1) ...
      / curTotalSpinesInSurround;
    
    
    % NOTE(amotta): Contrary to the above analysis, we now take into
    % account the observed size-dependence of spine density on spine size.
    curDensCtrl = curCellData.densCtrl;
    assert(abs(sum(curDensCtrl(:)) - 1) < 1E-3);
    
    curDensDiff = curDensReal - curDensCtrl;
    curDensDiffs(curDendIdx, :, :) = curDensDiff;
    
    curCondReal = curDensReal;
    curCondReal(~curCondReal) = eps;
    curCondReal = curCondReal ./ sum(curCondReal, 1);
    
    curCondCtrl = curDensCtrl;
    curCondCtrl(~curCondCtrl) = eps;
    curCondCtrl = curCondCtrl ./ sum(curCondCtrl, 1);
    
    curCondDiff = curCondReal - curCondCtrl;
    curCondDiffs(curDendIdx, :, :) = curCondDiff;
    curCondDiffDiag = curCondDiff(1:(imSize + 1):end);
    curCondDiffDiags(curDendIdx, :) = curCondDiffDiag;
end

curX = linspace(lims(1), lims(2), imSize);
curTitle = 'Grand average across proximal dendrites';


% Spine density plot
curMean = sum(curSpineCounts(:) .* curSpineDensDiffs, 1);
curMean = curMean / sum(curSpineCounts);
% NOTE(amotta): For plotting only.
curMask = curSpineCounts >= 100;

curFig = figure();
curAx = axes(curFig);
axis(curAx, 'square');
hold(curAx, 'on');

plot(curAx, ...
    curX, curSpineDensDiffs(curMask, :), ...
    'Color', repelem(0.7, 3));
plot(curAx, curX, zeros(size(curX)), 'k--');
plot(curAx, curX, curMean, 'k', 'LineWidth', 2);

xlabel(curAx, sprintf('log10(%s)', sizeLabel));
ylabel(curAx, 'ΔP(SH in surround) to random');
ylim(curAx, [-1, +1] * max(abs(ylim(curAx))));
if ~isempty(curTitle); title(curAx, curTitle); end
connectEM.Figure.config(curFig, info);
curFig.Position(3:4) = [440, 420];


% Size clustering plot
curMean = sum(curWeights(:) .* curCondDiffDiags, 1);
curMean = curMean / sum(curWeights);
% NOTE(amotta): For plotting only.
curMask = curWeights >= 100;

curFig = figure();
curAx = axes(curFig);
axis(curAx, 'square');
hold(curAx, 'on');

plot(curAx, ...
    curX, curCondDiffDiags(curMask, :), ...
    'Color', repelem(0.7, 3));
plot(curAx, curX, zeros(size(curX)), 'k--');
plot(curAx, curX, curMean, 'k', 'LineWidth', 2);

xlabel(curAx, sprintf('log10(%s)', sizeLabel));
ylabel(curAx, 'ΔP(SH of identical size in surround) to random');
curYlim = curCondDiffDiags(curMask, 75:225);
curYlim = max(abs(curYlim(:)));
ylim(curAx, [-1, +1] * curYlim);

if ~isempty(curTitle); title(curAx, curTitle); end
connectEM.Figure.config(curFig, info);
curFig.Position(3:4) = [440, 420];


% Joint size relationship
curMean = sum(curWeights(:) .* curDensDiffs, 1);
curMean = shiftdim(curMean, 1) / sum(curWeights);

curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');

imagesc(curAx, curMean);
curCbar = colorbar(curAx);

xlabel(curAx, curXLabel);
ylabel(curAx, curYLabel);
curAx.CLim = [-1, +1] * max(abs(curMean(:)));
curCbar.Label.String = 'ΔP(ref. SH, surround SH) to random';

curConfigAxis(curAx);
colormap(curAx, curBlueRedCmap);
plot(curAx, xlim(curAx), ylim(curAx), 'k--');
if ~isempty(curTitle); title(curAx, curTitle); end
connectEM.Figure.config(curFig, info);


% Conditional size relationship
curMean = sum(curWeights(:) .* curCondDiffs, 1);
curMean = shiftdim(curMean, 1) / sum(curWeights);

curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');

imagesc(curAx, curMean);
curCbar = colorbar(curAx);

xlabel(curAx, curXLabel);
ylabel(curAx, curYLabel);
curAx.CLim = [-1, +1] * max(abs(curMean(:)));
curCbar.Label.String = 'ΔP(surround SH | ref. SH) to random';

curConfigAxis(curAx);
colormap(curAx, curBlueRedCmap);
plot(curAx, xlim(curAx), ylim(curAx), 'k--');
if ~isempty(curTitle); title(curAx, curTitle); end
connectEM.Figure.config(curFig, info);


%% Utilities
function [im, bwOut] = kde2d(x, y, w, imSize, xlim, ylim, bwIn)
    bwOut = bwIn;
    
   [gridY, gridX] = ndgrid( ...
        linspace(ylim(1), ylim(2), imSize), ...
        linspace(xlim(1), xlim(2), imSize));
    grid = horzcat(gridY(:), gridX(:));

    optKvPairs = cell(0, 2);
    if ~isempty(bwIn); optKvPairs(end + 1, :) = {'Bandwidth', bwIn}; end
    if ~isempty(w); optKvPairs(end + 1, :) = {'Weights', w}; end
    optKvPairs = transpose(optKvPairs);
    
    im = cat(2, y(:), x(:));
    if isempty(im); im = nan(imSize); return; end
   [im, ~, bw] = ksdensity(im, grid, optKvPairs{:});
    im = reshape(im, [imSize, imSize]);
    
    if isempty(bwOut); bwOut = bw; end
end
