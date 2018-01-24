% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_18_a.mat');
outDir = '/tmpscratch/amotta/l4/2018-01-24-axons-18a-isosurfaces';

% temporarily add Benedikt's repo to path
beneDir = '/gaba/u/amotta/code/benedikt';
oldPath = addpath(genpath(beneDir), '-begin');
restorePath = onCleanup(@() path(oldPath));
clear oldPath;

info = Util.runInfo();

%% load parameters
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

% for Benedikt's code
param.agglo = struct;
param.agglo.axonAggloFile = axonFile;

% for speed
param.seg = struct;
param.seg.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
param.seg.backend = 'wkwrap';

%% load and complete agglomerates
cacheFile = fullfile(outDir, 'axons.mat');

if ~exist(cacheFile, 'file')
    disp('Generating axon agglomerates and cache');
    axons = L4.Axons.getLargeAxons(param, true, true);
    save(cacheFile, 'info', 'axons');
else
    disp('Loading axon agglomerates from cache');
    load(cacheFile, 'axons');
end

%% generate isosurfaces
Visualization.exportAggloToAmira( ...
    param, axons, outDir, 'reduce', 0.05, ...
    'smoothSizeHalf', 4, 'smoothWidth', 8);