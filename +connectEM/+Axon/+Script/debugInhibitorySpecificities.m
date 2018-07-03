% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
outputDir = '/home/amotta/Desktop/inh-sd-specific';

minSynPre = 10;

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);

[conn, syn, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

%% Preparing data
[conn, axonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, axonClasses, 'minSynPre', minSynPre);
axonClasses = ...
    connectEM.Connectome.buildAxonSpecificityClasses(conn, axonClasses);

%% Export
clear cur*;
mkdir(outputDir);

curTargetClass = 'SmoothDendrite';

rng(0);
axonIds = axonClasses(2).specs.(curTargetClass).axonIds;
axonIds = axonIds(randperm(numel(axonIds)));
axonIds = axonIds(1:min(25, numel(axonIds)));

synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.type = syn.synapses.type(synT.id);

numDigits = ceil(log10(1 + numel(axonIds)));

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

for curIdx = 1:numel(axonIds)
    curId = conn.axonMeta.id(axonIds(curIdx));
    curSynT = synT(synT.preAggloId == curId, :);
    curSynT.targetClass = conn.denMeta.targetClass(curSynT.postAggloId);
    
    curAxon = conn.axons(curId);
    curSynapses = cellfun( ...
        @vertcat, ...
        syn.synapses.presynId(curSynT.id), ...
        syn.synapses.postsynId(curSynT.id), ...
        'UniformOutput', false);
    
    curNodes = [curAxon; curSynapses];
    curNodes = cellfun( ...
        @(segIds) segPoints(segIds, :), ...
        curNodes, 'UniformOutput', false);
    
    curSkel = Skeleton.fromMST( ...
        curNodes, param.raw.voxelSize, skel);
    
    curSkel.names{1} = sprintf('Axon %d', curId);
    curSkel.colors{1} = [0, 0, 1, 1];
    
    curSkel.names(2:end) = arrayfun( ...
        @(id, type, target) sprintf( ...
            'Synapse %d (%s onto %s)', id, type, target), ...
        curSynT.id, curSynT.type, curSynT.targetClass, ...
        'UniformOutput', false);
    curSkel.colors(2:end) = {[1, 1, 0, 1]};
    
    curSkelName = sprintf( ...
        '%0*d_axon-%d.nml', numDigits, curIdx, curId);
    curSkel.write(fullfile(outputDir, curSkelName));
end
