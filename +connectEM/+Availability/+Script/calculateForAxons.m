% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
blockDataFile = '/tmpscratch/amotta/l4/2018-02-02-surface-availability-connectome-axons-18-a/block-data.mat';
axonsPerJob = 10;

outFile = fileparts(blockDataFile);
outFile = fullfile(outFile, 'axon-availability.mat');

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

blockData = load(blockDataFile);

%% capture target availability for every µm
maxDist = 1 + diff(param.bbox, 1, 2);
maxDist = maxDist .* param.raw.voxelSize(:);
maxDist = ceil(sqrt(sum(maxDist .^ 2)));

saveDists = 0:1E3:maxDist;

%% calculate on cluster
axonCount = numel(blockData.axonBlocks);

inputArgs = arrayfun(@(i) {{ ...
    i:min(i + axonsPerJob - 1, axonCount)}}, ...
    1:axonsPerJob:axonCount);
sharedInputArgs = {param, blockDataFile, saveDists};

cluster = Cluster.getCluster( ...
    '-pe openmp 1', ...
    '-l h_vmem=12G', ...
    '-l h_rt=12:00:00', ...
    '-l s_rt=11:59:00');

job = Cluster.startJob( ...
    @jobFunction, inputArgs, ...
    'sharedInputs', sharedInputArgs, ...
    'numOutputs', 1, ...
    'cluster', cluster, ....
    'name', mfilename);
wait(job);

%% build output
classCount = numel(blockData.targetClasses);

outMat = matfile(outFile);
outMat.axonAvail = nan(classCount, numel(saveDists), axonCount);

tic;
fprintf('Writing output file... ');
for curTaskIdx = 1:numel(inputArgs)
    curTask = job.Tasks(curTaskIdx);
    
    curAxonIds = inputArgs{curTaskIdx}{1};
    curResults = curTask.OutputArguments{1};
    
    % write partial results to output file
    outMat.axonAvail(:, :, curAxonIds) = curResults;
end

fprintf('done!');
toc;

%% runs on cluster
function targetClassAvail = jobFunction( ...
        param, blockDataFile, saveDists, axonIds)
    blockData = load(blockDataFile);
    
    targetClassAvail = ...
        connectEM.Availability.calculateForAxons( ...
            param, blockData, saveDists, axonIds);
end