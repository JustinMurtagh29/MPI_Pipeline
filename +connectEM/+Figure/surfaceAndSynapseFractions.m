% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
availBlockFile = '/tmpscratch/amotta/l4/2018-07-18-surface-availability-for-connectome-v7-partially-split/block-data.mat';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

boxCount = 100;
boxSizeVx = [480, 480, 240];

targetClasses = { ...
    'Somata', 'Soma'; ...
    'ProximalDendrite', 'PD'; ...
    'SmoothDendrite', 'SD'; ...
    'ApicalDendrite', 'AD'; ...
    'AxonInitialSegment', 'AIS'; ...
    'Rest', 'Rest'};

targetLabels = targetClasses(:, 2);
targetClasses = targetClasses(:, 1);

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);

[conn, syn] = connectEM.Connectome.load(param, connFile);
conn = connectEM.Connectome.prepareForSpecificityAnalysis(conn);

availBlock = load(availBlockFile);
availBlockSize = availBlock.info.param.blockSize;
assert(~any(mod(boxSizeVx, availBlockSize)));

%% Prepare availabilities
[~, sortIds] = ismember( ...
    availBlock.targetClasses, ...
    targetClasses(1:(end - 1)));
sortIds(~sortIds) = ...
    numel(targetClasses):numel(sortIds);

% Build 'Rest' class
availBlock = availBlock.targetClassAreas;
availBlock(sortIds, :, :, :) = availBlock;
availBlock(numel(targetClasses), :, :, :) = ...
    sum(availBlock(numel(targetClasses):end, :, :, :), 1);
availBlock = availBlock(1:numel(targetClasses), :, :, :);

%% Prepare synapses
synT = connectEM.Connectome.buildSynapseTable(conn, syn);

% Determine target class id
synT.targetClassId = ...
    conn.denMeta.targetClass(synT.postAggloId);
[~, synT.targetClassId] = ismember( ...
    synT.targetClassId, targetClasses(1:(end - 1)));
synT.targetClassId(~synT.targetClassId) = numel(targetClasses);
assert(all(synT.targetClassId));

% Calculate position
synT.pos = cellfun( ...
    @vertcat, ...
    syn.synapses.presynId(synT.id), ...
    syn.synapses.postsynId(synT.id), ...
    'UniformOutput', false);
synT.pos = cellfun( ...
    @(segIds) mean(segPoints(segIds, :), 1), ...
    synT.pos, 'UniformOutput', false);
synT.pos = cell2mat(synT.pos);

%% Sample random boxes
segBoxSize = 1 + diff(param.bbox, 1, 2)';

rng(0);
randBoxes = arrayfun( ...
    @(n) randi(n, boxCount, 1), ...
    floor(segBoxSize ./ boxSizeVx), ...
    'UniformOutput', false);
randBoxes = cell2mat(randBoxes);
randBoxes = (randBoxes - 1) .* boxSizeVx;
randBoxes = param.bbox(:, 1)' + randBoxes;

%% Analyse boxes
boxData = nan(numel(targetClasses), 2, boxCount);

for curIdx = 1:boxCount
    curBox = randBoxes(curIdx, :);
    curBox = [curBox; curBox + boxSizeVx]; %#ok
    
    % Surface fraction
    curAvailIds = curBox - param.bbox(:, 1)';
    curAvailIds = curAvailIds ./ availBlockSize + 1;
    
    curSurfFrac = availBlock(:, ...
            curAvailIds(1, 1):curAvailIds(2, 1), ...
            curAvailIds(1, 2):curAvailIds(2, 2), ...
            curAvailIds(1, 3):curAvailIds(2, 3));
	curSurfFrac = reshape(curSurfFrac, size(curSurfFrac, 1), []);
    curSurfFrac = sum(curSurfFrac, 2) ./ sum(curSurfFrac(:));
    
    % Synapse fraction
    curSynT = synT( ...
        all(synT.pos >= curBox(1, :), 2) ...
      & all(synT.pos  < curBox(2, :), 2), :);
    
    curSynFrac = accumarray( ...
        curSynT.targetClassId, 1, size(targetClasses));
    curSynFrac = curSynFrac ./ sum(curSynFrac);
    
    boxData(:, 1, curIdx) = curSurfFrac;
    boxData(:, 2, curIdx) = curSynFrac;
end

% Get rid of boxes that don't have any synapses
dropMask = squeeze(any(any(isnan(boxData), 1), 2));
boxData(:, :, dropMask) = [];

%% Global means
globalSurfFrac = sum(reshape(availBlock, numel(targetClasses), []), 2);
globalSurfFrac = globalSurfFrac ./ sum(globalSurfFrac);

globalSynFrac = accumarray( ...
    synT.targetClassId, 1, size(targetClasses));
globalSynFrac = globalSynFrac ./ sum(globalSynFrac);

globalData = cat(2, globalSurfFrac, globalSynFrac);

%% Show results
clear cur*;

curDataX = reshape(boxData(:, 1, :), [], 1);
curDataY = reshape(boxData(:, 2, :), [], 1);

curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [540, 375];

curAx = axes(curFig);
axis(curAx, 'square');
hold(curAx, 'on');

curAx.XLim = [0.001, 1];
curAx.YLim = [0.001, 1];

for curIdx = 1:numel(targetClasses)
    curColor = curAx.ColorOrder(curIdx, :);
    curDataX = reshape(boxData(curIdx, 1, :), [], 1);
    curDataY = reshape(boxData(curIdx, 2, :), [], 1);
    
    scatter(curAx, curDataX, curDataY, 3 * 36, curColor, '.');
end

for curIdx = 1:numel(targetClasses)
    curColor = curAx.ColorOrder(curIdx, :);
    curPlot = scatter(curAx, ...
        globalData(curIdx, 1), ...
        globalData(curIdx, 2), ...
        2 * 36, curColor, ...
        'o', 'filled');
    curPlot.MarkerEdgeColor = 'black';
end

curAx.Box = 'off';
curAx.TickDir = 'out';
curAx.XScale = 'log';
curAx.YScale = 'log';

xlabel(curAx, 'Surface fraction');
ylabel(curAx, 'Synapse fraction');

curPlots = flip(cat(1, curAx.Children));

[curLeg, curIcons] = legend( ...
    curPlots(1:numel(targetClasses)), ...
    targetClasses, 'Location', 'eastoutside');
curLeg.Box = 'off';

curMarkers = curIcons((end - numel(targetClasses) + 1):end);
curMarkers = cat(1, curMarkers.Children);
set(curMarkers, 'MarkerSize', 24);

line(curAx, curAx.XLim, curAx.YLim, 'Color', 'black');

annotation( ...
    curFig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
	'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');
