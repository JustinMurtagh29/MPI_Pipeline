% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered_classified.mat');
outputDir = '/home/amotta/Desktop/spine-shaft-axons';

minSynCount = 10;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);

[conn, syn] = ...
    connectEM.Connectome.load( ...
        param, connFile, synFile);

%% Complete axon meta data
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.type = syn.synapses.type(synT.id);

axonMeta = conn.axonMeta;
axonMeta.priSpineSynCount = accumarray( ...
    synT.preAggloId, synT.type == 'PrimarySpine', size(conn.axons));
axonMeta.priSpineSynFrac = axonMeta.priSpineSynCount ./ axonMeta.synCount;

axonMeta.shaftSynCount = accumarray( ...
    synT.preAggloId, synT.type == 'Shaft', size(conn.axons));
axonMeta.shaftSynFrac = axonMeta.shaftSynCount ./ axonMeta.synCount;

% Get rid of axons with too few synapses
axonMeta(axonMeta.synCount < minSynCount, :) = [];

%% Export examples
rng(0);
spineAxonIds = find(axonMeta.priSpineSynFrac >= 0.7);
spineAxonIds = spineAxonIds(randperm(numel(spineAxonIds)));
spineAxonIds = reshape(spineAxonIds(1:25), [], 1);

rng(0);
shaftAxonIds = find(axonMeta.shaftSynFrac >= 0.7);
shaftAxonIds = shaftAxonIds(randperm(numel(shaftAxonIds)));
shaftAxonIds = reshape(shaftAxonIds(1:25), [], 1);

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

axonIds = {spineAxonIds; shaftAxonIds};
axonTag = repelem({'spine'; 'shaft'}, cellfun(@numel, axonIds));
axonIds = cell2mat(axonIds);

numDigits = ceil(log10(1 + numel(axonIds)));

for curIdx = 1:numel(axonIds)
    curId = axonMeta.id(axonIds(curIdx));
    curSynT = synT(synT.preAggloId == curId, :);
    curTag = axonTag{curIdx};
    
    curAxon = conn.axons(curId);
    curSynapses = syn.synapses.postsynId(curSynT.id);
    
    curNodes = [curAxon; curSynapses];
    curNodes = cellfun( ...
        @(segIds) segPoints(segIds, :), ...
        curNodes, 'UniformOutput', false);
    
    curSkel = Skeleton.fromMST( ...
        curNodes, param.raw.voxelSize, skel);
    
    curSkel.names{1} = sprintf('Axon %d', curId);
    curSkel.colors{1} = [0, 0, 1, 1];
    
    curSkel.names(2:end) = arrayfun( ...
        @(id) sprintf('Synapse %d', id), ...
        curSynT.id, 'UniformOutput', false);
    curSkel.colors(2:end) = {[1, 1, 0, 1]};
    
    curSkelName = sprintf( ...
        '%0*d_%s_axon-%d.nml', ...
        numDigits, curIdx, curTag, curId);
    curSkel.write(fullfile(outputDir, curSkelName));
end
