% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered_classified.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_spine_attachment.mat');
outputDir = '/home/amotta/Desktop/thalamocortical';

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);

[conn, syn] = connectEM.Connectome.load(param, connFile, synFile);

shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

%% Prepare synapse list
shLUT = Agglo.buildLUT(maxSegId, shAgglos);
synT = connectEM.Connectome.buildSynapseTable(conn, syn);

synT.shId = cellfun( ...
    @(segIds) max(shLUT(segIds)), ...
    syn.synapses.postsynId(synT.id));
synT(~synT.shId, :) = [];

%% Export to webKNOSSOS
rng(0);
axonIds = find(conn.axonMeta.axonClass == 'Thalamocortical');
axonIds = axonIds(randperm(numel(axonIds)));
axonIds = axonIds(1:20);

numDigits = ceil(log10(1 + numel(axonIds)));

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

for curIdx = 1:numel(axonIds)
    curId = axonIds(curIdx);
    curAxon = conn.axons(curId);
    
    curSynT = synT(synT.preAggloId == curId, :);
    curSyns = shAgglos(curSynT.shId);
    
    curNodes = cellfun( ...
        @(segIds) segPoints(segIds, :), ...
        cat(1, curAxon, curSyns), ...
        'UniformOutput', false);
    curSkel = Skeleton.fromMST( ...
        curNodes, param.raw.voxelSize, skel);
    
    curSkel.names{1} = sprintf('Axon %d', curId);
    curSkel.colors{1} = [0, 0, 1, 1];
    
    curSkel.names(2:end) = arrayfun( ...
        @(id) sprintf('Spine head %d', id), ...
        curSynT.shId, 'UniformOutput', false);
    curSkel.colors(2:end) = {[1, 1, 0, 1]};
    
    curSkelName = sprintf( ...
        '%0*d_axon-%d.nml', numDigits, curIdx, curId);
    curSkel.write(fullfile(outputDir, curSkelName));
end

%% Export axons only for MH (but all of them)
axonIds = find(conn.axonMeta.axonClass == 'Thalamocortical');
numDigits = ceil(log10(1 + numel(axonIds)));

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

nodes = cellfun( ...
	@(segIds) segPoints(segIds, :), ...
	conn.axons(axonIds), 'UniformOutput', false);

skel = Skeleton.fromMST( ...
    nodes, param.raw.voxelSize, skel);
skel.names = arrayfun( ...
    @(id) sprintf('Axon %d', id), ...
    axonIds, 'UniformOutput', false);

skel.write(fullfile(outputDir, 'all.nml'));

