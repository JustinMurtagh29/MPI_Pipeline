% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a.mat');
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3.mat');

outputDir = '/home/amotta/Desktop/ad-specific-axons';
assert(mkdir(outputDir));

info = Util.runInfo();

%% loading data
syn = load(synFile);
conn = load(connFile);

param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

points = Seg.Global.getSegToPointMap(param);

%% build class connectome
% try to replicate cass connectome
connectome = conn.connectome;

% count synapses
connectome.synCount = cellfun(@numel, connectome.synIdx);
connectome.synIdx = [];

% add target class to connectome
[~, denMetaRow] = ismember(connectome.edges(:, 2), conn.denMeta.id);
connectome.targetClass = conn.denMeta.targetClass(denMetaRow);

[targetClasses, ~, connectome.targetClassId] = ...
    unique(connectome.targetClass);

classConnectome = accumarray( ...
    cat(2, connectome.edges(:, 1), connectome.targetClassId), ...
    connectome.synCount, [numel(conn.axons), numel(targetClasses)]);

%% building skeletons
className = 'ApicalDendrite';

conn.axonMeta.spineSynFrac = ...
    conn.axonMeta.spineSynCount ...
    ./ conn.axonMeta.synCount;

axonSpecs = ...
    classConnectome ...
    ./ sum(classConnectome, 2);

axonIds = find( ...
    conn.axonMeta.synCount >= 10 ...
  & abs(conn.axonMeta.spineSynFrac - 0.5) > 0.2 ...
  & axonSpecs(:, targetClasses == className) > 0.45);

% sort by increasing spine fraction
[~, sortIds] = sort(conn.axonMeta.spineSynFrac(axonIds));
axonIds = axonIds(sortIds);

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', mfilename, info.git_repos{1}.hash));

% targets
dendIds = find(conn.denMeta.targetClass == className);
for curIdx = 1:numel(dendIds)
    curDendId = dendIds(curIdx);
    curAgglo = conn.dendrites{curDendId};
    
    skel = Skeleton.fromMST( ...
        points(curAgglo, :), param.raw.voxelSize, skel);
    skel.names{end} = sprintf('Dendrite %d', curDendId);
    skel.colors{end} = [0, 0, 1, 1];
end

% axons
for curIdx = 1:numel(axonIds)
    curAxonId = axonIds(curIdx);
    curAgglo = conn.axons{curAxonId};
    
    % axon
    curSkel = Skeleton.fromMST( ...
        points(curAgglo, :), param.raw.voxelSize, skel);
    curSkel.names{end} = sprintf( ...
        'Axon %d. Spine fraction %.2f', ...
        curAxonId, conn.axonMeta.spineSynFrac(curAxonId));
    curSkel.colors{end} = [1, 0, 0, 1];
    
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
    
    curSkel = Skeleton.fromMST( ...
        curSyns, param.raw.voxelSize, curSkel);
    
    curSynTreeIds = (1:numel(curSynIds)) - 1;
    curSynTreeIds = numel(curSkel.nodes) - curSynTreeIds;
    curSkel.names(curSynTreeIds) = arrayfun( ...
        @(id) sprintf('Synapse %d', id), ...
        curSynIds, 'UniformOutput', false);
    curSkel.colors(curSynTreeIds) = {[1, 1, 0, 1]};
    
    curSkelFile = sprintf('%0*d_axon-%d.nml', ...
        ceil(log10(1 + numel(axonIds))), curIdx, curAxonId);
    curSkel.write(fullfile(outputDir, curSkelFile));
end