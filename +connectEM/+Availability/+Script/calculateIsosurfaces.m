% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outDir = '/tmpscratch/amotta/l4/2018-04-30-availabilities-soma-axon/iso';

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

distVol = struct;
distVol.rootDir = '/tmpscratch/amotta/l4/2018-04-30-availabilities-soma-axon/wkw-lz4';
distVol.thresh = 5000;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

% for speed
param.seg = struct;
param.seg.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
param.seg.backend = 'wkwrap';

conn = load(connFile);

%% Select dendrite subset
postSynCount = accumarray( ...
    conn.connectome.edges(:, 2), ...
    cellfun(@numel, conn.connectome.synIdx)', ...
   [numel(conn.dendrites), 1]);

postSynIds = find(postSynCount >= minSynCount);
postSynAgglos = conn.dendrites(postSynIds);

%% Building isosurface
Visualization.exportAggloToAmira( ...
    param, postSynAgglos, outDir, ...
    'distVol', distVol, 'reduce', 0.05, ...
    'smoothSizeHalf', 4, 'smoothWidth', 8);
