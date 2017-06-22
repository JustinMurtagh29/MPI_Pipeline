% Parameters
outputFolder = '/gaba/scratch/mberning/axonQueryGenerationReverse/';

% Load data
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
[graph, segmentMeta, borderMeta, globalSegmentPCA] = connectEM.loadAllSegmentationData(p);

% Save state for bookkeeping
load(['/gaba/scratch/mberning/axonQueryGeneration/beforeQueryGeneration.mat'], 'directionality', 'axonsNew');

% Generate axon queries
querySaveFolder = [outputFolder 'queriesMat/'];
if ~exist(querySaveFolder, 'dir')
    mkdir(querySaveFolder)
end

options.latentScore = 0.7;
options.segDirScore = 0.9;
options.border = [3000; -3000];
options.writeTasksToFile = true;
options.boundingBoxForTasks = false;
options.reverse = true;

connectEM.generateAxonQueries(p, graph, segmentMeta, borderMeta, directionality, axonsNew, querySaveFolder, options);

