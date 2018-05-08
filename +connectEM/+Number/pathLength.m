% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
trunkFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3.mat');

[~, outFile] = fileparts(connFile);
outFile = sprintf('%s_pathLengths.mat', outFile);
outFile = fullfile(fileparts(connFile), outFile);

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = connectEM.Connectome.load(param, connFile);

dendrites = load(conn.info.param.dendriteFile);
dendrites = dendrites.dendrites;

trunks = load(trunkFile);
trunks = trunks.dendrites(trunks.indBigDends);

voxelSize = param.raw.voxelSize;
maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);

%% Axons
Util.log('Calculating axon path lengths');

axonPathLengths = SuperAgglo.fromAgglo( ...
    conn.axons, segPoints, 'mst', 'voxelSize', voxelSize);
axonPathLengths = SuperAgglo.pathLength(axonPathLengths, voxelSize);

%% Dendrites
Util.log('Calculating dendrite path lengths');

dendritePathLengths = SuperAgglo.pathLength( ...
    dendrites(conn.denMeta.parentId), voxelSize);

%% Dendrite trunks (prior to spine attachment)
Util.log('Calculating dendritic trunk path lengths');

trunkIds = Agglo.buildLUT(maxSegId, Agglo.fromSuperAgglo(trunks));
trunkIds = cellfun(@(ids) mode(nonzeros(trunkIds(ids))), conn.dendrites);
trunkMask = ~isnan(trunkIds);

trunkPathLengths = nan(size(trunkIds));
trunkPathLengths(trunkMask) = ...
    SuperAgglo.pathLength(trunks(trunkIds(trunkMask)), voxelSize);

%% Writing results
out = struct;
out.axonPathLengths = axonPathLengths;
out.dendritePathLengths = dendritePathLengths;
out.trunkPathLengths = trunkPathLengths;
out.trunkIds = trunkIds;
out.info = info;

Util.saveStruct(outFile, out);
Util.protect(outFile);
