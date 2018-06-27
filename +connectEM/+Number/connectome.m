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

[conn, syn] = connectEM.Connectome.load(param, connFile);

%% Numbers
numAxons = numel(conn.axons) %#ok
numDendrites = numel(conn.dendrites) %#ok

axonMask = conn.axonMeta.synCount >= 10;
numAxonsWithTenSynapses = sum(axonMask) %#ok

dendMask = conn.denMeta.synCount >= 10;
numDendritesWithTenSynapses = sum(dendMask) %#ok

numSomata = sum(...
    conn.denMeta.targetClass(dendMask) == 'Somata') %#ok
numSmoothDendrites = sum( ...
    conn.denMeta.targetClass(dendMask) == 'SmoothDendrite') %#ok
numApicalDendrites = sum( ...
    conn.denMeta.targetClass(dendMask) == 'ApicalDendrite') %#ok
numAxonInitialSegments = sum( ...
    conn.denMeta.targetClass(dendMask) == 'AxonInitialSegment') %#ok
numAxonInitialSegmentsAll = sum( ...
    conn.denMeta.targetClass == 'AxonInitialSegment') %#ok
