% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
outDir = '/home/amotta/Desktop/inh-axons-seeded-at-ad-shaft-syns';

minSynPre = 10;

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

points = Seg.Global.getSegToPointMap(param);

[conn, syn, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

%% Prepare data
% Copy & paste from connectEM.Specificity.plot
[conn, axonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, axonClasses, 'minSynPre', minSynPre);
[classConn, targetClasses] = ...
    connectEM.Connectome.buildClassConnectome(conn);

% Build list of possible seed synapses
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.ontoTargetClass = conn.denMeta.targetClass(synT.postAggloId);
synT.type = syn.synapses.type(synT.id);

%% Select seed synapses
inhAxonIds = axonClasses(2).axonIds;

seedSynT = synT( ...
    ismember(synT.preAggloId, inhAxonIds) ...
  & synT.ontoTargetClass == 'ApicalDendrite' ...
  & synT.type == 'Shaft', :);

% NOTE(amotta): Instead of randomly sampling, let's just shuffle the entire
% table of seed synapses and then "select" the first N rows.
rng(0);
seedSynT = seedSynT(randperm(height(seedSynT)), :);

%% Export to NML files
% Copy & Paste from connectEM.Connectome.debugInhibitoryAxonsAndSynapses
% with minor modifification for this specific purpose.
exportIds = 1:50;

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

mkdir(outDir);

for curIdx = exportIds
    curAxonId = seedSynT.preAggloId(curIdx);
    curAxonSegIds = conn.axons{curAxonId};
    
    curSeedSynId = seedSynT.id(curIdx);
    curSeedSynSegIds = unique(vertcat( ...
        syn.synapses.presynId{curSeedSynId}, ...
        syn.synapses.postsynId{curSeedSynId}));
    
    curSkel = skel;
    curSkel = Skeleton.fromMST( ...
        points(curAxonSegIds, :), param.raw.voxelSize, curSkel);
    curSkel.names{end} = sprintf('Axon %d', curAxonId);
    curSkel.colors{end} = [0, 0, 1, 1];
    
    curSkel = Skeleton.fromMST( ...
        points(curSeedSynSegIds, :), ...
        param.raw.voxelSize, curSkel);
    curSkel.names{end} = sprintf( ...
        'Synapse %d (seed)', curSeedSynId);
    curSkel.colors{end} = [1, 0, 0, 1];
    
    curSkelName = sprintf( ...
        '%0*d_axon-%d.nml', ...
        ceil(log10(1 + numel(exportIds))), ...
        curIdx, curAxonId);
    curSkel.write(fullfile(outDir, curSkelName));
end
