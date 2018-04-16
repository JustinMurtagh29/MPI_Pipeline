% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
aggloFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_v5.mat');
outDir = '/tmpscratch/amotta/l4/2018-04-16-axon-initial-segment-isosurfaces';

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

agglos = load(aggloFile);
agglos = agglos.dendrites(agglos.indAIS);
agglos = Agglo.fromSuperAgglo(agglos);

% for speed
param.seg = struct;
param.seg.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
param.seg.backend = 'wkwrap';

%% Generate isosurfaces
Visualization.exportAggloToAmira( ...
    param, agglos, outDir, 'reduce', 0.05, ...
    'smoothSizeHalf', 4, 'smoothWidth', 8);
