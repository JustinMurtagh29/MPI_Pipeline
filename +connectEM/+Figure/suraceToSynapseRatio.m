% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
availBlockFile = '/tmpscratch/amotta/l4/2018-02-02-surface-availability-connectome-axons-18-a/block-data.mat';
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');

boxCount = 10;
boxSizeVx = [480, 480, 240];

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);

syn = load(synFile);
conn = load(connFile);

availBlock = load(availBlockFile);
availBlockSize = availBlock.info.param.blockSize;
assert(~any(mod(boxSizeVx, availBlockSize)));

targetClasses = ...
    unique(conn.denMeta.targetClass);
[~, targetClassIds] = ismember( ...
    targetClasses, availBlock.targetClasses);
assert(all(targetClassIds));

availBlock = ...
    availBlock.targetClassAreas( ...
        targetClassIds, :, :, :);

%% Prepare synapses
synT = connectEM.Connectome.buildSynapseTable(conn, syn);

% Determine target class id
synT.targetClassId = conn.denMeta.targetClass(synT.postAggloId);
[~, synT.targetClassId] = ismember(synT.targetClassId, targetClasses);

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
for curIdx = 1:boxCount
    curBox = randBoxes(curIdx, :);
    curBox = [curBox; curBox + boxSizeVx]; %#ok
    
    curSynT = synT( ...
        all(synT.pos >= curBox(1, :), 2) ...
      & all(synT.pos  < curBox(2, :), 2), :);
    curSynFrac = accumarray( ...
        curSynT.targetClassId, 1, size(targetClasses));
    curSynFrac = curSynFrac ./ sum(curSynFrac);
    
    curAvailIds = curBox - param.bbox(:, 1)';
    curAvailIds = curAvailIds ./ availBlockSize + 1;
    
    curAvailFrac = availBlock(:, ...
            curAvailIds(1, 1):curAvailIds(2, 1), ...
            curAvailIds(1, 2):curAvailIds(2, 2), ...
            curAvailIds(1, 3):curAvailIds(2, 3));
	curAvailFrac = reshape(curAvailFrac, size(curAvailFrac, 1), []);
    curAvailFrac = sum(curAvailFrac, 2) ./ sum(curAvailFrac(:));
end
