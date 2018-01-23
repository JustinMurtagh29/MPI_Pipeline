% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_16_b.mat');
outDir = '/tmpscratch/amotta/l4/2018-02-23-all-axon-isosurfaces';

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

%% load and complete agglomerates
axons = L4.Axons.getLargeAxons(param, true, true);
save(fullfile(outDir, 'axons.mat'), 'info', 'axons');
