% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
aggloFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');
outDir = '/tmpscratch/amotta/l4/2018-04-16-distance-volume-test/20180417-isosurfaces';

distVol = struct;
distVol.rootDir = '/tmpscratch/amotta/l4/2018-04-16-distance-volume-test/wkw';
distVol.thresh = 10000;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

% for speed
param.seg = struct;
param.seg.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
param.seg.backend = 'wkwrap';

agglos = load(aggloFile);
agglos = agglos.dendrites;

%% Render single soma for prototyping
% TODO(amotta): Remove this section once it's working.
maxSegId = Seg.Global.getMaxSegId(param);
aggloLUT = Agglo.buildLUT(maxSegId, agglos);

aggloId = Seg.Global.getSegIds(param, [1950, 6423, 858]);
aggloId = aggloLUT(aggloId);
agglos = agglos(aggloId);

%% Building isosurface
Visualization.exportAggloToAmira( ...
    param, agglos, outDir, ...
    'distVol', distVol, 'reduce', 0.05, ...
    'smoothSizeHalf', 4, 'smoothWidth', 8);
