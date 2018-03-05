% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');
debugDir = '/home/amotta/Desktop/inhibitory-axons';

info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

points = Seg.Global.getSegToPointMap(param);

syn = load(synFile);
conn = load(connFile);

%% find inhibitory axons
conn.axonMeta.spineSynFrac = ...
    conn.axonMeta.spineSynCount ...
    ./ conn.axonMeta.synCount;

axonIds = find( ...
    conn.axonMeta.synCount >= 10 ...
  & conn.axonMeta.spineSynFrac <= 0.3);

%% export examples to webKNOSSOS
rng(0);
axonIds = axonIds(randperm(numel(axonIds)));

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

for curIdx = 1
    curAxonId = axonIds(curIdx);
    curAxonSegIds = conn.axons{curAxonId};
    
    curSynIds = (conn.connectome.edges(:, 1) == curAxonId);
    curSynIds = cell2mat(conn.connectome.synIdx(curSynIds));
    
    curSynSegIds = cellfun( ...
        @(a, b) unique(vertcat(a, b)), ...
        syn.synapses.presynId(curSynIds), ...
        syn.synapses.postsynId(curSynIds), ...
        'UniformOutput', false);
    
    curSkel = skel;
    curSkel = Skeleton.fromMST( ...
        points(curAxonSegIds, :), param.raw.voxelSize, curSkel);
    curSkel.names{end} = sprintf('Axon %d', curAxonId);
    curSkel.colors{end} = [0, 0, 1, 1];
    
    for curSynIdx = 1:numel(curSynIds)
        curSkel = Skeleton.fromMST( ...
            points(curSynSegIds{curSynIdx}, :), ...
            param.raw.voxelSize, curSkel);
        curSkel.names{end} = sprintf( ...
            'Synapse %d', curSynIds(curSynIdx));
        curSkel.colors{end} = [1, 0, 0, 1];
    end
end