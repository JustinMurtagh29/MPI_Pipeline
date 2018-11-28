% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
outputDir = '/home/amotta/Desktop/exc-specs';

minSynPre = 10;
info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);

[conn, syn, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

[conn, axonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, axonClasses, 'minSynPre', minSynPre);

%% Prepare data
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.ontoTarget = conn.denMeta.targetClass(synT.postAggloId);

specClasses = ...
    connectEM.Connectome.buildAxonSpecificityClasses( ...
        conn, axonClasses(1:2));

%% Export samples
clear cur*;

rng(0);
mkdir(outputDir);

curAxonIds = specClasses(1).specs.ProximalDendrite.axonIds;
curAxonIds = curAxonIds(randperm(numel(curAxonIds)));
curAxonIds = curAxonIds(1:10);

curTemp = skeleton();
curTemp = Skeleton.setParams4Pipeline(curTemp, param);
curTemp = curTemp.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

curNumDigits = ceil(log10(1 + numel(curAxonIds)));

for curIdx = 1:numel(curAxonIds)
    curId = conn.axonMeta.id(curAxonIds(curIdx));
    curSynT = synT(synT.preAggloId == curId, :);
    curSynT.type = syn.synapses.type(curSynT.id);
    
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
        curNodes, param.raw.voxelSize, curTemp);
    
    curSkel.names{1} = sprintf('Axon %d', curId);
    curSkel.colors{1} = [0, 0, 1, 1];
    
    curSkel.names(2:end) = arrayfun( ...
        @(id, type, target) sprintf( ...
            'Synapse %d (%s onto %s)', id, type, target), ...
        curSynT.id, curSynT.type, curSynT.ontoTarget, ...
        'UniformOutput', false);
    curSkel.colors(2:end) = {[1, 1, 0, 1]};
    
    curSkelName = sprintf( ...
        '%0*d_axon-%d.nml', curNumDigits, curIdx, curId);
    curSkel.write(fullfile(outputDir, curSkelName));
end
