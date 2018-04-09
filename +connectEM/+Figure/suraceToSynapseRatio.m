% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
availBlockFile = '/tmpscratch/amotta/l4/2018-02-02-surface-availability-connectome-axons-18-a/block-data.mat';
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');

availBlockFile = '/home/amotta/Desktop/cache/block-data.mat';
synFile = '/home/amotta/Desktop/cache/SynapseAgglos_v3_ax_spine_clustered.mat';
connFile = '/home/amotta/Desktop/cache/connectome_axons_18_a_ax_spine_syn_clust.mat';

boxCount = 1;
boxSizeVx = 15 * 32;
boxSizeVx = [boxSizeVx, boxSizeVx, boxSizeVx / 2];

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);

availBlock = load(availBlockFile);
availBlockSize = availBlock.info.param.blockSize;
assert(~any(mod(boxSizeVx, availBlockSize)));

syn = load(synFile);
conn = load(connFile);

targetClasses = ...
    unique(conn.denMeta.targetClass);
[~, targetClassIds] = ismember( ...
    targetClasses, availBlock.targetClasses);
assert(all(targetClassIds));

%% Calculate synapse location
synapses = syn.synapses;
synapses.pos = cellfun( ...
    @vertcat, ...
    synapses.presynId, ...
    synapses.postsynId, ...
    'UniformOutput', false);
synapses.pos = cell2mat(cellfun( ...
    @(segIds) mean(segPoints(segIds, :), 1), ...
    synapses.pos, 'UniformOutput', false));

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

%%
for curIdx = 1:boxCount
    curBox = randBoxes(curIdx, :);
    curBox = [curBox; curBox + boxSizeVx];
    
    curSynIds = find( ...
        all(synapses.pos >= curBox(1, :), 2) ...
      & all(synapses.pos  < curBox(2, :), 2));
    
    curAvailIds = curBox - param.bbox(:, 1)';
    curAvailIds = curAvailIds ./ availBlockSize + 1;
    
    curAvail = ...
        availBlock.targetClassAreas(:, ...
            curAvailIds(1, 1):curAvailIds(2, 1), ...
            curAvailIds(1, 2):curAvailIds(2, 2), ...
            curAvailIds(1, 3):curAvailIds(2, 3));
    
	curAvail = reshape(curAvail, size(curAvail, 1), []);
    curAvail = sum(curAvail, 2) ./ sum(curAvail(:));
end
