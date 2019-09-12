% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connOldFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');
connNewFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

adAxonIdOld = 23325;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

connOld = load(connOldFile);
[connNew, synNew] = connectEM.Connectome.load(param, connNewFile);
synNewT = connectEM.Connectome.buildSynapseTable(connNew, synNew);

axonsNew = load(connNew.info.param.axonFile);
axonsNew = axonsNew.axons(axonsNew.indBigAxons);

maxSegId = Seg.Global.getMaxSegId(param);
segPos = Seg.Global.getSegToPointMap(param);

%% Translate indices to old connectome
axonNewLUT = Agglo.buildLUT(maxSegId, connNew.axons);

adAxonIdNew = cellfun( ...
    @(ids) mode(nonzeros(axonNewLUT(ids))), ...
    connOld.axons(adAxonIdOld));

%% Evaluate synapses
[classConn, targetClasses] = ...
    connectEM.Connectome.buildClassConnectome(connNew);

disp(targetClasses);
disp(classConn(adAxonIdNew, :));

%% Build NML of axon with synapses
clear cur*;

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = Skeleton.setDescriptionFromRunInfo(skel, info);

curAxon = connNew.axonMeta.parentId(adAxonIdNew);
curAxon = axonsNew(curAxon);

curSynT = synNewT;
curSynT = synNewT(curSynT.preAggloId == adAxonIdNew, :);
curSynT.targetClass = connNew.denMeta.targetClass(curSynT.postAggloId);

curSynPos = synNew.synapses;
curSynPos = curSynPos(curSynT.id, :);

curSynPos = cellfun( ...
    @vertcat, ...
    curSynPos.presynId, ...
    curSynPos.postsynId, ...
    'UniformOutput', false);
curSynPos = cellfun( ...
    @(ids) segPos(ids, :), curSynPos, ...
    'UniformOutput', false);

skel = skel.addTree( ...
    sprintf('Axon %d', adAxonIdNew), ...
    curAxon.nodes(:, 1:3), curAxon.edges);
skel = Skeleton.fromMST( ...
    curSynPos, param.raw.voxelSize, skel);
skel.names(2:end) = arrayfun( ...
    @char, curSynT.targetClass, ...
    'UniformOutput', false);
