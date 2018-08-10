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
[connNew, synNew, axonClasses] = ...
    connectEM.Connectome.load(param, connNewFile);

maxSegId = Seg.Global.getMaxSegId(param);

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
