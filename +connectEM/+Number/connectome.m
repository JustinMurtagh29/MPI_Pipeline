% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, axonClasses] = connectEM.Connectome.load(param, connFile);

%% Connectome
numAxons = numel(conn.axons) %#ok
numDendrites = numel(conn.dendrites) %#ok

axonMask = conn.axonMeta.synCount >= 10;
numAxonsWithTenSynapses = sum(axonMask) %#ok

dendMask = conn.denMeta.synCount >= 10;
aisMask = conn.denMeta.targetClass == 'AxonInitialSegment';
numDendritesWithTenSynapses = sum(dendMask) %#ok

numSomata = sum(...
    conn.denMeta.targetClass(dendMask) == 'Somata') %#ok
numSmoothDendrites = sum( ...
    conn.denMeta.targetClass(dendMask) == 'SmoothDendrite') %#ok
numApicalDendrites = sum( ...
    conn.denMeta.targetClass(dendMask) == 'ApicalDendrite') %#ok
numAxonInitialSegments = sum(dendMask & aisMask) %#ok
numAxonInitialSegmentsAll = sum(aisMask) %#ok

%% Axon classes
numLikelyExcitatoryAxons = numel(axonClasses(1).axonIds) %#ok
numLikelyInhibitoryAxons = numel(axonClasses(2).axonIds) %#ok

fractionOfAxonsLikelyBeingExcitatory = ...
    numLikelyExcitatoryAxons / ( ...
    numLikelyExcitatoryAxons + numLikelyInhibitoryAxons) %#ok

%% Synapse fractions
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.targetClass = conn.denMeta.targetClass(synT.postAggloId);

synT.inConnectome = ...
    axonMask(synT.preAggloId) & ( ...
    dendMask(synT.postAggloId) | ...
    aisMask(synT.postAggloId));

numSynpases = height(synT) %#ok
numSomaSynapses = sum(synT.targetClass == 'Somata') %#ok
fractionOfSynapsesOntoSomata = mean(synT.targetClass == 'Somata') %#ok

numSynapsesInConnectome = sum(synT.inConnectome) %#ok
numSomaSynapsesInConnectome = sum( ...
    synT.targetClass(synT.inConnectome) == 'Somata') %#ok
fractionOfSynapsesInCommectomeOntoSomata = mean( ...
    synT.targetClass(synT.inConnectome) == 'Somata') %#ok

% Analyse likely inhibitory synapses
synT.isLikelyInhibitory = ismember( ...
    synT.preAggloId, axonClasses(2).axonIds);
numLikelyInhSynapses = sum(synT.isLikelyInhibitory) %#ok
numLikelyInhSomaSynapses = sum( ...
    synT.targetClass(synT.isLikelyInhibitory) == 'Somata') %#ok
fractionOfLikelyInhibitorySynapsesOntoSomata = mean( ...
    synT.targetClass(synT.isLikelyInhibitory) == 'Somata') %#ok
