% Search for shaft-preferring excitatory axons.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
outputDir = '/home/amotta/Desktop';

minSynPre = 10;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);

[conn, syn, axonClasses] = connectEM.Connectome.load(param, connFile);
[synIds, synsInConn] = connectEM.Axon.getSynapses(conn, syn);

%% Prepare data
conn.axonMeta.fullPriSpineSynFrac = ...
    conn.axonMeta.fullPriSpineSynCount ...
 ./ conn.axonMeta.fullSynCount;

debugConfigs = struct;
debugConfigs(1).axonIds = axonClasses(3).axonIds;
debugConfigs(1).tag = 'thalamocortical';

%% Search for potential shaft-preferring excitatory axons
clear cur*;

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

inConnStr = {', not in conn.', ''};

for curConfig = debugConfigs
    curAxonIds = curConfig.axonIds;
    curOutputDir = fullfile(outputDir, curConfig.tag);
    
    assert(~exist(curOutputDir, 'dir'));
    mkdir(curOutputDir);

    rng(0);
    curAxonIds = curAxonIds(randperm(numel(curAxonIds)));
    curAxonIds = curAxonIds(1:50);

    curNumDigits = ceil(log10(1 + numel(curAxonIds)));

    for curIdx = 1:numel(curAxonIds)
        curId = conn.axonMeta.id(curAxonIds(curIdx));
        
        curSynIds = synIds{curId};
        curSynsInConn = synsInConn{curId};
        curSynTypes = syn.synapses.type(curSynIds);
        
        curAxon = conn.axons(curId);
        curSynapses = cellfun( ...
            @vertcat, ...
            syn.synapses.presynId(curSynIds), ...
            syn.synapses.postsynId(curSynIds), ...
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
                'Synapse %d (%s%s)', ...
                id, type, inConnStr{1 + inConn}), ...
            curSynIds, curSynTypes, curSynsInConn, ...
            'UniformOutput', false);
        curSkel.colors(2:end) = {[1, 1, 0, 1]};

        curSkelName = sprintf( ...
            '%0*d_axon-%d.nml', curNumDigits, curIdx, curId);
        curSkel.write(fullfile(curOutputDir, curSkelName));
    end
end
