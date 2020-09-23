% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

outputDir = '/home/amotta/Desktop/ad-specific-axons';
assert(mkdir(outputDir));

info = Util.runInfo();
Util.showRunInfo(info);

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

points = Seg.Global.getSegToPointMap(param);

[conn, syn, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);
[conn, axonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, axonClasses(1:2), 'minSynPre', 10);
axonClasses = ...
    connectEM.Connectome.buildAxonSpecificityClasses(conn, axonClasses);

%% building skeletons
clear cur*;

axonIds = axonClasses(2).specs.ApicalDendrite.axonIds;
dendIds = find(conn.denMeta.targetClass == 'ApicalDendrite');

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = Skeleton.setDescriptionFromRunInfo(skel, info);

% targets
for curIdx = 1:numel(dendIds)
    curDendId = dendIds(curIdx);
    curAgglo = conn.dendrites{curDendId};
    
    skel = Skeleton.fromMST( ...
        points(curAgglo, :), param.raw.voxelSize, skel);
    skel.names{end} = sprintf('Dendrite %d', curDendId);
    skel.colors{end} = [0, 0, 1, 1];
end

[skel, curGroupId] = skel.addGroup('Dendrites');
skel = skel.addTreesToGroup(1:numel(dendIds), curGroupId);

% axons
for curIdx = 1:numel(axonIds)
    curAxonId = axonIds(curIdx);
    curAgglo = conn.axons{curAxonId};
    
    % axon
    skel = Skeleton.fromMST( ...
        points(curAgglo, :), param.raw.voxelSize, skel);
    skel.names{end} = sprintf('Axon %d', curAxonId);
    skel.colors{end} = [1, 0, 0, 1];
    
    % synapses
    curSynIds = (conn.connectome.edges(:, 1) == curAxonId);
    curSynIds = cell2mat(conn.connectome.synIdx(curSynIds));
    
    curSyns = cellfun( ...
        @vertcat, ...
        syn.synapses.presynId(curSynIds), ...
        syn.synapses.postsynId(curSynIds), ...
        'UniformOutput', false);
    curSyns = cellfun( ...
        @(segIds) points(unique(segIds), :), ...
        curSyns, 'UniformOutput', false);
    
    skel = Skeleton.fromMST( ...
        curSyns, param.raw.voxelSize, skel);
    
    curSynTreeIds = (1:numel(curSynIds)) - 1;
    curSynTreeIds = numel(skel.nodes) - curSynTreeIds;
    skel.names(curSynTreeIds) = arrayfun( ...
        @(id) sprintf('Synapse %d', id), ...
        curSynIds, 'UniformOutput', false);
    skel.colors(curSynTreeIds) = {[1, 1, 0, 1]};
    
    % group axon and synapses
    curTreeIds = 1:(numel(curSynIds) + 1);
    curTreeIds = numel(skel.nodes) + 1 - curTreeIds;
   [skel, curGroupId] = skel.addGroup(sprintf('Axon %d', curAxonId));
    skel = skel.addTreesToGroup(curTreeIds, curGroupId);
end

skel.write(fullfile(outputDir, 'skel.nml'));
