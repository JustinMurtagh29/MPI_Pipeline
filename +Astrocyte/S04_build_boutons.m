% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

[connDir, connName] = fileparts(connFile);

interSynFile = sprintf('%s_intersynapse_v2.mat', connName);
interSynFile = fullfile(connDir, interSynFile);

outFile = sprintf('%s_axonalBoutons_v1.mat', connName);
outFile = fullfile(connDir, outFile);
clear connDir connName;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
% param.seg = struct;
% param.seg.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
% param.seg.backend = 'wkwrap';

segPoints = Seg.Global.getSegToPointMap(param);
segPoints = param.raw.voxelSize .* segPoints;

segVols = Seg.Global.getSegToSizeMap(param);
segVols = prod(param.raw.voxelSize) .* segVols;

conn = load(connFile);
syn = load(conn.info.param.synFile);
interSyn = load(interSynFile);

%% Build bouton agglomerates and calculate volume
boutonAgglos = connectEM.Axon.buildBoutonAgglos( ...
    segPoints, conn, syn, interSyn, 'parallelize', true);

boutonVols = cellfun( ...
    @(agglos) cellfun( ...
        @(segIds) sum(segVols(segIds)), agglos), ...
    boutonAgglos, 'UniformOutput', false);

%% Build output
out = struct;
out.info = info;
out.boutonAgglos = boutonAgglos;
out.boutonVols = boutonVols;

Util.saveStruct(save('/gaba/u/yyener/astrocyte/synapses/boutons.mat', out);
Util.protect(outFile);
