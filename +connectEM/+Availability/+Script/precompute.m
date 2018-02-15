% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
blockSize = [32, 32, 16];
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a.mat');
outDir = '/tmpscratch/amotta/l4/2018-02-09-surface-availability-connectome-axons-18-a';

info = Util.runInfo();

%% load parameter
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

%% precompute
[axonBlocks, targetClassAreas, targetClasses] = ...
	connectEM.Availability.precompute(param, connFile, blockSize);

%% store result
outFile = fullfile(outDir, 'block-data.mat');
Util.save(outFile, info, axonBlocks, targetClassAreas, targetClasses);
