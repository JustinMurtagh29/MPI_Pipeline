% Search for shaft-preferring excitatory axons.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
outputDir = '/home/amotta/Desktop';

minSynPre = 10;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);

[conn, syn] = connectEM.Connectome.load(param, connFile);

synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.type = syn.synapses.type(synT.id);

%% Prepare data
conn.axonMeta.fullPriSpineSynFrac = ...
    conn.axonMeta.fullPriSpineSynCount ...
 ./ conn.axonMeta.fullSynCount;

axonClasses = struct;
axonClasses(1).axonIds = find( ...
    conn.axonMeta.synCount >= minSynPre ...
  & conn.axonMeta.fullPriSpineSynFrac <= 0.2);
axonClasses(1).tag = 'less-than-20-percent-pri-spines';

axonClasses(2).axonIds = find( ...
    conn.axonMeta.synCount >= minSynPre ...
  & conn.axonMeta.fullPriSpineSynFrac > 0.2 ...
  & conn.axonMeta.fullPriSpineSynFrac < 0.5);
axonClasses(2).tag = 'between-20-and-50-percent-pri-spines';

%% Search for potential shaft-preferring excitatory axons
clear cur*;

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

for curAxonClass = axonClasses
    curAxonIds = curAxonClass.axonIds;
    curOutputDir = fullfile(outputDir, curAxonClass.tag);
    
    assert(~exist(curOutputDir, 'dir'));
    mkdir(curOutputDir);

    rng(0);
    curAxonIds = curAxonIds(randperm(numel(curAxonIds)));
    curAxonIds = curAxonIds(1:25);

    curNumDigits = ceil(log10(1 + numel(curAxonIds)));

    for curIdx = 1:numel(curAxonIds)
        curId = conn.axonMeta.id(curAxonIds(curIdx));
        curSynT = synT(synT.preAggloId == curId, :);

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
            @(id, type) sprintf('Synapse %d (%s)', id, type), ...
            curSynT.id, curSynT.type, 'UniformOutput', false);
        curSkel.colors(2:end) = {[1, 1, 0, 1]};

        curSkelName = sprintf( ...
            '%0*d_axon-%d.nml', curNumDigits, curIdx, curId);
        curSkel.write(fullfile(curOutputDir, curSkelName));
    end
end
