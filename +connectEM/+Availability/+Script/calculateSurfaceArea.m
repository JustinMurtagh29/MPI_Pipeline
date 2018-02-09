% Written by
%   Kevin M. Boergens <kevin.boergens@brain.mpg.de>
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
blockSize = [32, 32, 16];
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a.mat');
outDir = '/tmpscratch/amotta/l4/2018-02-02-surface-availability-connectome-axons-18-a';

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

%% actually start job
inputArgs = arrayfun(@(i) {{i}}, 1:numel(param.local));
sharedInputArgs = {param, blockSize, connFile, outDir};
mkdir(outDir);

cluster = Cluster.getCluster( ...
    '-pe openmp 1', ...
    '-l h_vmem=24G', ...
    '-l h_rt=1:00:00', ...
    '-l s_rt=0:59:00');

job = Cluster.startJob( ...
    @connectEM.Availability.calculateSurfaceArea, ...
    inputArgs, 'sharedInputs', sharedInputArgs, ...
    'cluster', cluster, 'name', mfilename);
