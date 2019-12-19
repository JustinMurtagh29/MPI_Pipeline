% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
trunkFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3.mat');
dendFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto-and-manual.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
asiRunId = '20190227T082543';

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

%% Let's have a look at synapses onto proximal dendrites
clear cur*;
rng(0);

% NOTE(amotta): The first run corresponds to actually observed data. All
% subsequent runs use randomly redistributed ASI areas.
curNumRuns = 1 + 100;
curDistThresh = 2000;

curShT = table;
curAsiT = asiT( ...
    asiT.type == 'PrimarySpine' ...
  & asiT.targetClass == 'ProximalDendrite' ...
  & ismember(asiT.axonClass, {'Corticocortical', 'Thalamocortical'}), :);
curAsiT.postAggloId = shT.dendId(curAsiT.shId);

seedAreas = nan([height(curAsiT), curNumRuns]);
condAreas = cell([height(curAsiT), curNumRuns]);

for curRunIdx = 1:curNumRuns    
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

        curOtherShMask = curDendDists < curDistThresh;
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

%% Plot results
curDens = cell(curNumRuns, 1);
curLims = [-1.5, 0];

curConfigAxis = @(ax) set(ax, ...
    'PlotBoxAspectRatio', [1, 1, 1], ...
    'DataAspectRatioMode', 'auto', ...
    'XDir', 'normal', 'YDir', 'normal');

for curRunIdx = 1:curNumRuns
    curSeedAreas = seedAreas(:, curRunIdx);
    curCondAreas = condAreas(:, curRunIdx);
    
    curX = cellfun(@numel, curCondAreas);
    curX = repelem(curSeedAreas, curX);
    curY = cell2mat(curCondAreas);

    curX = log10(curX);
    curY = log10(curY);

   [~, curDens{curRunIdx}] = ...
        connectEM.Libs.kde2d( ...
            cat(2, curX, curY), 256, ...
            repelem(curLims(1), 2), ...
            repelem(curLims(2), 2));
    curDens{curRunIdx} = curDens{curRunIdx} / sum(curDens{curRunIdx}(:));
end

curRealDens = curDens{1};
curCtrlDens = mean(cat(3, curDens{2:end}), 3);
curDiffDens = curRealDens - curCtrlDens;

curFig = figure();
curAx = axes(curFig);
imagesc(curAx, curRealDens);
curConfigAxis(curAx);
connectEM.Figure.config(curFig, info);

curFig = figure();
curAx = axes(curFig);
imagesc(curAx, curDiffDens);
curAx.CLim = [-1, +1] * max(abs(curDiffDens(:)));
curConfigAxis(curAx);
connectEM.Figure.config(curFig, info);

curRealCond = sum(curRealDens, 1);
curRealCond(~curRealCond) = eps;
curRealCond = curRealDens ./ curRealCond;

curCtrlCond = sum(curCtrlDens, 1);
curCtrlCond(~curCtrlCond) = eps;
curCtrlCond = curCtrlDens ./ curCtrlCond;
curDiffCond = curRealCond - curCtrlCond;

curFig = figure();
curAx = axes(curFig);
imagesc(curAx, curRealCond);
curConfigAxis(curAx);
connectEM.Figure.config(curFig, info);

curFig = figure();
curAx = axes(curFig);
imagesc(curAx, curDiffCond);
curClim = curDiffCond(:, 50:end);
curClim = [-1, +1] * max(abs(curClim(:)));
curAx.CLim = curClim;
curConfigAxis(curAx);
connectEM.Figure.config(curFig, info);
