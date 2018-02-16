% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
shFile = fullfile(rootDir, 'aggloState', 'spine_heads_and_attachment_03.mat');
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');

% set this variable to debug
debugDir = '/home/amotta/Desktop/multi-axon-sh';

info = Util.runInfo();

%% loading data
conn = load(connFile);
syn = load(synFile);

shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

% for debugging
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);
points = Seg.Global.getSegToPointMap(param);

%% assign synapses to spine heads
shLUT = Agglo.buildLUT(maxSegId, shAgglos);

% limit to spine synapses
synT = syn.synapses;
synT.id = reshape(1:size(synT, 1), [], 1);
synT.shId = cellfun(@(ids) max(shLUT(ids)), syn.synapses.postsynId);

axonT = table;
[axonT.id, ~, uniRows] = unique(conn.connectome.edges(:, 1));
axonT.synIds = accumarray( ...
    uniRows, reshape(1:size(conn.connectome, 1), [], 1), ...
    [], @(rows) {cell2mat(conn.connectome.synIdx(rows))});
axonT.shIds = cellfun(@(ids) ...
    reshape(setdiff(synT.shId(ids), 0), [], 1), ...
    axonT.synIds, 'UniformOutput', false);

shIds = cell2mat(axonT.shIds);
axonIds = repelem(axonT.id, cellfun(@numel, axonT.shIds));

shT = table;
[shT.id, ~, uniRows] = unique(shIds);
shT.axonIds = accumarray(uniRows, axonIds, [], @(ids) {ids(:)});
shT(cellfun(@isscalar, shT.axonIds), :) = [];

%% export random examples
rng(0);
randIds = randperm(size(shT, 1), 10);
randShT = shT(randIds, :);

for curIdx = 1:size(randShT, 1)
    curSh = randShT(curIdx, :);
    curSh = table2struct(curSh);
    
    curNodes = [ ...
        shAgglos(curSh.id);
        conn.axons(curSh.axonIds)];
    curNodes = cellfun( ...
        @(segIds) points(segIds, :), ...
        curNodes, 'UniformOutput', false);
    
    curNames = arrayfun( ...
        @(axonId) sprintf('Axon #%d', axonId), ...
        curSh.axonIds, 'UniformOutput', false);
    curNames = [ ...
        {sprintf('Spinehead #%d', curSh.id)}; curNames];
    
    curSkel = Skeleton.fromMST( ...
        curNodes, param.raw.voxelSize);
    curSkel.names = curNames;
    
    curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    
    curSkelFile = sprintf( ...
        '%d_spinehead-%d.nml', curIdx, curSh.id);
    curSkel.write(fullfile(debugDir, curSkelFile));
end