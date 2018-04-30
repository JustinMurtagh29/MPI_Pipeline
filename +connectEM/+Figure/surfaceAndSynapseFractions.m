% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
availBlockFile = '/tmpscratch/amotta/l4/2018-04-27-surface-availability-connectome-v5/block-data.mat';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

boxCount = 100;
boxSizeVx = [480, 480, 240];

targetClasses = { ...
    'Somata', 'Soma'; ...
    'ApicalDendrite', 'AD'; ...
    'SmoothDendrite', 'SD'; ...
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
boxGroups = reshape(1:numel(targetClasses), [], 1);
boxGroups = repmat(boxGroups, 1, 2, size(boxData, 3));

boxGroups(:, 1, :) = 2 .* boxGroups(:, 1, :) - 1;
boxGroups(:, 2, :) = 2 .* boxGroups(:, 2, :);

% Add jitter
boxJitter = 0.5 * (rand(size(boxGroups)) - 0.5);
boxJitter = boxGroups + boxJitter;

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [720, 300];

ax = axes(fig);
hold(ax, 'on');

colors = ax.ColorOrder([1, 3], :);
colorGroup = repmat(1:2, 1, numel(targetClasses));

% Fake plots for legend
plot(ax, nan, nan, '.', 'Color', colors(1, :));
plot(ax, nan, nan, '.', 'Color', colors(2, :));

scatter(ax, ...
    boxJitter(:), boxData(:), [], ...
    colors(1 + mod(boxGroups(:) - 1, 2), :), '.');

for curIdx = 1:numel(globalData)
    curY = globalData(curIdx);
    curY = repelem(curY, 1, 2);
    
    curGroup = boxGroups(curIdx);
    curColor = colors(1 + mod(curGroup - 1, 2), :);
    curX = curGroup + [-0.45, +0.45];
    
    plot(ax, ...
        curX, curY, ...
        'Color', curColor, ...
        'LineWidth', 2);
end

ax.Box = 'off';
ax.TickDir = 'out';
ax.YScale = 'log';

ylim(ax, [1e-4, 1]);
yticklabels(ax, arrayfun( ...
    @(f) sprintf('%g', f), ...
    yticks(ax), 'UniformOutput', false));
ax.YAxis.MinorTick = 'off';
ylabel(ax, 'Fraction');

xticks(ax, 2 .* (1:numel(targetClasses)) - 0.5);
xticklabels(ax, targetLabels);

h = legend(ax, ...
    'Surface fraction', ...
    'Synapse fraction', ...
    'Location', 'eastoutside');
h.Box = 'off';

title(ax, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
