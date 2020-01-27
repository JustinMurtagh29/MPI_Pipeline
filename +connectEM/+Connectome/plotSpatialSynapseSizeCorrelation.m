% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
trunkFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3.mat');
dendFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto-and-manual.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
asiRunId = '20190227T082543';

cellDataFile = '/tmpscratch/amotta/l4/2020-01-27-cell-based-homeostatic-plasticity-analysis';
cellDataFile = fullfile(cellDataFile, 'cell-data_weighted_2-5um-radius_v1.mat');

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
imSize = 256;
lims = [-1.5, 0];
bw = [0.1, 0.1];

% NOTE(amotta): The first run corresponds to actually observed data. All
% subsequent runs use randomly redistributed ASI areas.
numRuns = 1 + 10;
distThresh = 2500;

if ~exist(cellDataFile, 'file')
    curAggloToCell = conn.denMeta.cellId;
    cellIds = setdiff(curAggloToCell, 0);
    cellData = cell(numel(cellIds), 1);

    parfor curCellIdx = 1:numel(cellIds)
        curCellId = cellIds(curCellIdx);

        curShT = shT;
        curAsiT = asiT( ...
            asiT.type == 'PrimarySpine' ...
          & asiT.targetClass == 'ProximalDendrite', :); %#ok

        curAsiT.postCellId = curAggloToCell(curAsiT.postAggloId); %#ok
        curAsiT = curAsiT(curAsiT.postCellId == curCellId, :);

        curCellSeedAreas = nan([height(curAsiT), numRuns]);
        curCellCondAreas = cell([height(curAsiT), numRuns]);

        rng(0);
        for curRunIdx = 1:numRuns
            curShT.asiArea(:) = nan;
            curShT.asiArea(curAsiT.shId) = curAsiT.area;

            for curIdx = 1:height(curAsiT)
                curSeedShId = curAsiT.shId(curIdx);
                curSeedDendId = shT.dendId(curSeedShId);

                curSeedAsiArea = curShT.asiArea(curSeedShId);

                % Find other spine heads onto same dendrite
                curDendShIds = dendT.shIds{curSeedDendId}; %#ok
                curSeedShMask = curDendShIds == curSeedShId;
                assert(sum(curSeedShMask) == 1);

                % Restrict to other spine heads in surround
                curDendDists = dendT.shToShDists{curSeedDendId};
                curDendDists = curDendDists(:, curSeedShMask);

                curOtherShMask = curDendDists < distThresh;
                curOtherShMask = curOtherShMask & not(curSeedShMask);

                curOtherShIds = curDendShIds(curOtherShMask);
                curOtherAsiAreas = curShT.asiArea(curOtherShIds);
                curOtherAsiAreas = rmmissing(curOtherAsiAreas);

                curCellSeedAreas(curIdx, curRunIdx) = curSeedAsiArea;
                curCellCondAreas{curIdx, curRunIdx} = curOtherAsiAreas(:);
            end

            % Shuffle areas
            curRandIds = randperm(height(curAsiT));
            curAsiT.area = curAsiT.area(curRandIds);
        end

        curDens = cell(numRuns, 1);
        for curRunIdx = 1:numRuns
            curSeedAreas = curCellSeedAreas(:, curRunIdx);
            curCondAreas = curCellCondAreas(:, curRunIdx);

            curN = cellfun(@numel, curCondAreas);
            curX = repelem(curSeedAreas, curN);
            curY = cell2mat(curCondAreas);

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
            curW = repelem(1 ./ curN, curN);

            curX = log10(curX);
            curY = log10(curY);

           [curDens{curRunIdx}, curRunBw] = kde2d( ...
                curX, curY, curW, imSize, lims, lims, bw);

            curDens{curRunIdx} = ...
                curDens{curRunIdx} ...
              / sum(curDens{curRunIdx}(:));
        end

        curRealDens = curDens{1};
        curCtrlDens = cat(3, curDens{2:end});
        curCtrlDens = mean(curCtrlDens, 3);

        curCellData = struct;
        curCellData.seedAreas = curCellSeedAreas(:, 1);
        curCellData.condAreas = curCellCondAreas(:, 1);
        curCellData.realDens = curRealDens;
        curCellData.ctrlDens = curCtrlDens;

        cellData{curCellIdx} = curCellData;
    end

    cellData = cat(1, cellData{:});

    % NOTE(amotta): To tell MATLAB that we don't use these after the `parfor`.
    clear curSeedAreas curCondAreas curRealDens;
    
    curOut = struct;
    curOut.info = info;
    curOut.cellIds = cellIds(:);
    curOut.cellData = cellData(:);
    
    Util.saveStruct(cellDataFile, curOut);
    Util.protect(cellDataFile);
else
    % NOTE(amotta): Load data from cache
   [cellIds, cellData] = Util.load(cellDataFile, 'cellIds', 'cellData');
end

%% Plot size versus number of neighboring spines
%{
clear cur*;

curX = log10(curSeedAreas(:, 1));
curY = curCondAreas(:, 1);
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

curTicks = [-1, 0];
curTickLabels = {'-1', '0'};
curTickIds = linspace(lims(1), lims(2), imSize);
[~, curTickIds] = arrayfun(@(t) min(abs(t - curTickIds)), curTicks);

curXLabel = 'log10(reference ASI area [µm²])';
curYLabel = 'log10(ASI area within %g µm [µm³])';
curYLabel = sprintf(curYLabel, distThresh / 1E3);

curRedCmap = connectEM.Figure.whiteRed(256);
curBlueRedCmap = connectEM.Figure.blueWhiteRed(256);

curConfigAxis = @(ax) set(ax, ...
    'PlotBoxAspectRatio', [1, 1, 1], ...
    'DataAspectRatioMode', 'auto', ...
    'XLim', [0, imSize] + 0.5, 'XDir', 'normal', ...
    'XTick', curTickIds, 'XTickLabel', curTickLabels, ...
    'YLim', [0, imSize] + 0.5, 'YDir', 'normal', ...
    'YTick', curTickIds, 'YTickLabel', curTickLabels);

curDiffCondDiags = nan(numel(cellData), imSize);

for curCellIdx = 1:numel(cellData)
    curCellId = cellIds(curCellIdx);
    curCellData = cellData(curCellIdx);
    
    curRealDens = curCellData.realDens;
    curCtrlDens = curCellData.ctrlDens;
    
    if any(isnan(curRealDens(:))); continue; end
    if any(isnan(curCtrlDens(:))); continue; end
    
    curRealCond = curRealDens;
    curRealCond(~curRealCond) = eps;
    curRealCond = curRealCond ./ sum(curRealCond, 1);
    
    curCtrlCond = curCtrlDens;
    curCtrlCond(~curCtrlCond) = eps;
    curCtrlCond = curCtrlCond ./ sum(curCtrlCond, 1);
    
    curDiffDens = curRealDens - curCtrlDens;
    curDiffCond = curRealCond - curCtrlCond;
    
    curDiffCondDiag = curDiffCond(1:(size(curDiffCond, 1) + 1):end);
    curDiffCondDiags(curCellIdx, :) = curDiffCondDiag;
    
    curTitle = sprintf('Synapses onto neuron %d', curCellId);

    % Plot 1
    %{
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
    %}
end

curDiffCondDiags(all(isnan(curDiffCondDiags), 2), :) = [];
figure; plot(linspace(0, 1.5, imSize), curDiffCondDiags, 'LineWidth', 2);
figure; plot(linspace(0, 1.5, imSize), mean(curDiffCondDiags, 1), 'k', 'LineWidth', 2);

%% Find size range with signs of homeostatic suppression in surround
clear cur*;
rng(0);

distThresh = 2500;
minSynCount = 200;

curVecX = -1.5:0.05:0;

curAsiT = asiT( ...
    asiT.type == 'PrimarySpine' ...
  & asiT.targetClass == 'ProximalDendrite' ...
  & ismember(asiT.axonClass, {'Corticocortical', 'Thalamocortical'}), :);

% IMPORTANT(amotta): Use `postAggloId` before transforming it!
curAsiT.postCellId = conn.denMeta.cellId(curAsiT.postAggloId);
curAsiT.postAggloId = shT.dendId(curAsiT.shId);

% NOTE(amotta): Restrict to cell with at least X synapses
[curCellIds, ~, curCellSynCounts] = unique(curAsiT.postCellId);
curCellSynCounts = accumarray(curCellSynCounts, 1);

curCellIds = curCellIds(curCellSynCounts >= minSynCount);
curCellIds = 32;
curMatY = nan(numel(curCellIds), numel(curVecX));

for curCellIdx = 1:numel(curCellIds)
    curCellId = curCellIds(curCellIdx);
    curCellAsiT = curAsiT(curAsiT.postCellId == curCellId, :);

    curShT = shT;
    curShT.asiArea(:) = nan;
    curShT.asiArea(curCellAsiT.shId) = curCellAsiT.area;

    curSeedAreas = nan([height(curCellAsiT), 1]);
    curCondAreas = cell([height(curCellAsiT), 1]);
    for curIdx = 1:height(curCellAsiT)
        curSeedShId = curCellAsiT.shId(curIdx);
        curSeedDendId = curCellAsiT.postAggloId(curIdx);

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

        curSeedAreas(curIdx) = curSeedAsiArea;
        curCondAreas{curIdx} = curOtherAsiAreas(:);
    end
    
    % To log10-space
    curSeedAreas = log10(curSeedAreas);
    curCondAreas = cellfun(@log10, curCondAreas, 'UniformOutput', false);

    curCondN = cellfun(@numel, curCondAreas);
    curCondAreas = cell2mat(curCondAreas);

    for curIdxX = 1:numel(curVecX)
        curX = curVecX(curIdxX);
        
        % NOTE(amotta): Uniform kernel
        curCondW = curSeedAreas - curX;
        curCondW = double(abs(curCondW) < 0.2);

        assert(isequal(size(curCondW), size(curCondN)));
        curCondW = repelem(curCondW ./ curCondN, curCondN);
        curCondW = curCondW / sum(curCondW(:));

        curY = sum(curCondAreas(:) .* curCondW(:));
        curY = curY - mean(curSeedAreas);
        
        curMatY(curCellIdx, curIdxX) = curY;
    end
end

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
%}

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
