% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');
outDir = '/tmpscratch/amotta/l4/2018-01-24-axons-18a-isosurfaces';

%% load parameters
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

axons = load(axonFile, 'axons');
axons = axons.axons;

% for speed
param.seg = struct;
param.seg.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
param.seg.backend = 'wkwrap';

%% generate isosurfaces
Visualization.exportAggloToAmira( ...
    param, axons, outDir, 'reduce', 0.05, ...
    'smoothSizeHalf', 4, 'smoothWidth', 8);
