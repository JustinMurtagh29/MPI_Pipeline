% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);
[conn, syn] = connectEM.Connectome.load(param, connFile);

%% Restrict to secondary spine innervations
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.type = syn.synapses.type(synT.id);
synT(synT.type ~= 'SecondarySpine', :) = [];

%% Export examples
rng(0);
curRandIds = randperm(height(synT));
curRandIds = curRandIds(1:25);

curSynT = synT(curRandIds, :);
curSynT.agglo = cellfun( ...
    @vertcat, ...
    syn.synapses.presynId(curSynT.id), ...
    syn.synapses.postsynId(curSynT.id), ...
    'UniformOutput', false);

curSkel = skeleton();
curSkel = Skeleton.setParams4Pipeline(curSkel, param);
curSkel = curSkel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

curAgglos = cellfun( ...
    @(segIds) segPoints(segIds, :), ...
    curSynT.agglo, 'UniformOutput', false);
curSkel = Skeleton.fromMST( ...
    curAgglos, param.raw.voxelSize, curSkel);
curSkel.names = arrayfun( ...
    @(id) sprintf('Synapse %d', id), ...
    curSynT.id, 'UniformOutput', false);
