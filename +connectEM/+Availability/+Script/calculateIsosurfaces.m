% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outDir = '/tmpscratch/amotta/l4/2018-04-30-availabilities-soma-axon';
isoDir = fullfile(outDir, 'iso');

% NOTE(amotta): This is an older connectome. It has the same dendrite
% agglomerates as `connectome_axons_18_a_ax_spine_syn_clust`, but different
% synapses. I'm reusing this version, because we want to render exactly
% the same subset of dendrites as shown in figure 1.
%
% For more details, see
% `+connectEM/+Connectome/exportPostsynToAmira.m`
% (commit hash 7f1fcc376674b03d7c4259d4cb2fe60a2ca5e636)
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a.mat');
minSynCount = 10;

% NOTE(amotta): For historic reasons the postsynaptic agglomerates are
% rendered from an other version of the connectome. But for target class
% identification we will use the latest results.
dendMetaFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

distVol = struct;
distVol.rootDir = '/tmpscratch/amotta/l4/2018-04-30-availabilities-soma-axon/wkw-lz4';
distVol.thresh = 5000;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);

% for speed
param.seg = struct;
param.seg.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
param.seg.backend = 'wkwrap';

conn = load(connFile);

% Load traget classes
dendMeta = connectEM.Connectome.load(param, dendMetaFile);
dendMeta.denMeta.agglo = dendMeta.dendrites;
dendMeta = dendMeta.denMeta;

%% Select dendrite subset
postSynCount = accumarray( ...
    conn.connectome.edges(:, 2), ...
    cellfun(@numel, conn.connectome.synIdx)', ...
   [numel(conn.dendrites), 1]);

postSynIds = find(postSynCount >= minSynCount);
postSynAgglos = conn.dendrites(postSynIds);

%% Build meta data
% Find target class for each postsynaptic agglomerate
targetClasses = Agglo.buildLUT(maxSegId, dendMeta.agglo);
targetClasses = cellfun(@(segIds) mode( ...
    nonzeros(targetClasses(segIds))), postSynAgglos);
targetClasses = dendMeta.targetClass(targetClasses);

[targetClasses, ~, targetClassFiles] = unique(targetClasses);
targetClasses = categories(targetClasses);

targetClassFiles = accumarray( ...
    targetClassFiles, reshape(1:numel(targetClassFiles), [], 1), [], ...
    @(aggloIds) {arrayfun(@(id) ...
        fullfile(isoDir, 'ply', sprintf('iso-%d.ply', id)), ...
        aggloIds, 'UniformOutput', false)});

%% Save meta data
out = struct;
out.info = info;
out.targetClasses = targetClasses;
out.targetClassFiles = targetClassFiles;

infoFile = fullfile(outDir, 'info.mat');
Util.saveStruct(infoFile, out);
Util.protect(infoFile);

%% Building isosurface
Visualization.exportAggloToAmira( ...
    param, postSynAgglos, isoDir, ...
    'distVol', distVol, 'reduce', 0.05, ...
    'smoothSizeHalf', 4, 'smoothWidth', 8);
