% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
blockSize = [32, 32, 16];
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
outDir = '/tmpscratch/amotta/l4/2018-07-18-surface-availability-for-connectome-v7-partially-split';

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
