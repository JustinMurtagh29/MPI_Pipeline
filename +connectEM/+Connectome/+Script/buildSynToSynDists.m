% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

info = Util.runInfo();

%% Complete configuration
[~, outFile] = fileparts(connFile);
outFile = sprintf('%s_synToSynDists.mat', outFile);
outFile = fullfile(fileparts(connFile), outFile);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn] = connectEM.Connectome.load(param, connFile);

% Dendrites in super-agglomerate format
dendrites = load(conn.info.param.dendriteFile);
dendrites = dendrites.dendrites(conn.denMeta.parentId);

% Segment meta data
segSizes = Seg.Global.getSegToSizeMap(param);
segPoints = Seg.Global.getSegToPointMap(param);

%% Calculate intersynapse distances for dendrites
Util.log('Calculating synapse distances along dendrites');

[dendSynToSynDists, dendSynIds] = ...
    connectEM.Connectome.buildSynToSynDist( ...
        conn, syn, dendrites, 'post', ...
        'voxelSize', param.raw.voxelSize, ...
        'segWeights', segSizes);

%% Convert axons into super-agglomerates
Util.log('Building axon super-agglomerates');

axons = SuperAgglo.fromAgglo( ...
    conn.axons, segPoints, 'mst', ...
    'voxelSize', param.raw.voxelSize);
SuperAgglo.check(axons);

%% Calculate intersynapse distances for axons
Util.log('Calculating synapse distances along axons');

[axonSynToSynDists, axonSynIds] = ...
    connectEM.Connectome.buildSynToSynDist( ...
        conn, syn, axons, 'pre', ...
        'voxelSize', param.raw.voxelSize, ...
        'segWeights', segSizes);

%% Build and save output
out = struct;
out.info = info;
out.axonSynToSynDists = axonSynToSynDists;
out.axonSynIds = axonSynIds;
out.dendSynToSynDists = dendSynToSynDists;
out.dendSynIds = dendSynIds;

Util.saveStruct(outFile, out);
Util.protect(outFile);
