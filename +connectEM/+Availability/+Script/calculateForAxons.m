% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
blockDataFile = '/tmpscratch/amotta/l4/2018-02-02-surface-availability-connectome-axons-18-a/proto-block-data.mat';
axonsPerJob = 10;

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

blockData = load(blockDataFile);

%% capture target availability for every µm
maxDist = 1 + diff(param.bbox, 1, 2);
maxDist = maxDist .* param.raw.voxelSize(:);
maxDist = ceil(sqrt(sum(maxDist .^ 2)));

saveDists = 1E3:1E3:maxDist;

%% calculate on cluster
axonCount = numel(blockData.axonBlocks);

inputArgs = arrayfun(@(i) {{ ...
    i:min(i + axonsPerJob - 1, axonCount)}}, ...
    1:axonsPerJob:axonCount);
sharedInputArgs = {param, blockDataFile, saveDists};

cluster = Cluster.getCluster( ...
    '-pe openmp 1', ...
    '-l h_vmem=12G', ...
    '-l h_rt=1:00:00', ...
    '-l s_rt=0:59:00');

job = Cluster.startJob( ...
    @jobFunction, inputArgs, ...
    'sharedInputs', sharedInputArgs, ...
    'numOutputs', 1, ...
    'cluster', cluster, ....
    'name', mfilename);
wait(job);

%% build output
out = fetchOutputs(job);
targetClassAvail = vertcat(out{:});

%% runs on cluster
function targetClassAvail = jobFunction( ...
        param, blockDataFile, saveDists, axonIds)
    blockData = load(blockDataFile);
    
    targetClassAvail = ...
        connectEM.Availability.calculateForAxons( ...
            param, blockData, saveDists, axonIds);
end