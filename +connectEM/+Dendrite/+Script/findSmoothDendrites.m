% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');

minSynCount = 10;
info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

points = Seg.Global.getSegToPointMap(param);

conn = load(connFile);
syn = load(synFile);

%% Based on spine synapse fraction
% This is justified according to Yunfeng's data
candIds = find( ...
    conn.denMeta.synCount >= minSynCount ...
  & conn.denMeta.targetClass ~= 'Somata');

spineSynFrac = ...
    conn.denMeta.spineSynCount ...
    ./ conn.denMeta.synCount;
candIds(spineSynFrac(candIds) > 0.5) = [];

%% Export to webKNOSSOS
rng(0);
randIds = randperm(numel(candIds));
randIds = reshape(candIds(randIds(1:25)), 1, []);

dendPoints = cellfun( ...
    @(segIds) points(segIds, :), ...
    conn.dendrites(randIds), ...
    'UniformOutput', false);

dendNames = arrayfun( ...
    @(idx, id) sprintf( ...
        '%0*d. Agglomerate %d', ...
        ceil(log10(1 + numel(randIds))), idx, id), ...
	1:numel(randIds), randIds, 'UniformOutput', false);

skel = Skeleton.fromMST(dendPoints, param.raw.voxelSize);
skel.names = dendNames;

skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));