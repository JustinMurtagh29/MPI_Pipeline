% Search for shaft-preferring excitatory axons.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
outputDir = '/home/amotta/Desktop';

minSynPre = 10;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);

[conn, syn] = connectEM.Connectome.load(param, connFile);

%% Find synapses not contained in connectome
axonLUT = Agglo.buildLUT(maxSegId, conn.axons);
syn.synapses.axonId = cellfun( ...
    @(ids) setdiff(axonLUT(ids), 0), ...
    syn.synapses.presynId, ...
    'UniformOutput', false);

syn.synapses.axonId(~cellfun(@isscalar, syn.synapses.axonId)) = {nan};
syn.synapses.axonId = cell2mat(syn.synapses.axonId);

syn.synapses.inConn = false(height(syn.synapses), 1);
syn.synapses.inConn(cell2mat(conn.connectome.synIdx)) = true;

%% Prepare data
conn.axonMeta.fullPriSpineSynFrac = ...
    conn.axonMeta.fullPriSpineSynCount ...
 ./ conn.axonMeta.fullSynCount;

axonClasses = struct;
axonClasses(1).axonIds = find( ...
    conn.axonMeta.synCount >= minSynPre ...
  & conn.axonMeta.fullPriSpineSynFrac > 0.2 ...
  & conn.axonMeta.fullPriSpineSynFrac < 0.5);
axonClasses(1).tag = 'between-20-and-50-percent-pri-spines';

%% Search for potential shaft-preferring excitatory axons
clear cur*;

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

inConnStr = {'. Not in connectome', ''};

for curAxonClass = axonClasses
    curAxonIds = curAxonClass.axonIds;
    curOutputDir = fullfile(outputDir, curAxonClass.tag);
    
    assert(~exist(curOutputDir, 'dir'));
    mkdir(curOutputDir);

    curNumDigits = ceil(log10(1 + numel(curAxonIds)));

    for curIdx = 1:numel(curAxonIds)
        curId = conn.axonMeta.id(curAxonIds(curIdx));
        curSynT = syn.synapses(syn.synapses.axonId == curId, :);

        curAxon = conn.axons(curId);
        curSynapses = cellfun( ...
            @vertcat, ...
            curSynT.presynId, ...
            curSynT.postsynId, ...
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
            @(id, type, inConn) sprintf( ...
                'Synapse %d (%s%s)', id, type, inConnStr{1 + inConn}), ...
            curSynT.id, curSynT.type, curSynT.inConn, ...
            'UniformOutput', false);
        curSkel.colors(2:end) = {[1, 1, 0, 1]};

        curSkelName = sprintf( ...
            '%0*d_axon-%d.nml', curNumDigits, curIdx, curId);
        curSkel.write(fullfile(curOutputDir, curSkelName));
    end
end
