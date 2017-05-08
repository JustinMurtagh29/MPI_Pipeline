
% load data
%load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
%segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat', 'voxelCount', 'centroid', 'box', 'maxSegId', 'cubeIdx', 'point');
temp = load('/gaba/scratch/kboerg/directCycle/cycle006.mat', 'axonsFinal', 'metrics');
agglos = temp.axonsFinal;
aggloSize = cellfun(@(x)sum(segmentMeta.voxelCount(x)), agglos);
[aggloSize, aggloIdx] = sort(aggloSize, 'descend');
agglos = agglos(aggloIdx);
agglosGtAxon2 = temp.axonsFinal(temp.metrics.axon1.foundAgglomerates_col{2});
clear temp;
% calculate queries and write to file on wk
outputFolder = '/gaba/scratch/mberning/flightQueries/';

% 2 different versions for MH (1 and 2 micron bounding box)
options.writeTasksToFile = true;
options.queryBoundingBoxSize = 1000;
options.datasetBorderExclusionSize = 2000;
q = connectEM.generateQueriesFromAgglos(p, segmentMeta, agglosGtAxon2, outputFolder, options);
connectEM.debugNewQueries(segmentMeta, agglosGtAxon2, q, outputFolder);
% 2nd version
options.writeTasksToFile = true;
options.queryBoundingBoxSize = 2000;
options.datasetBorderExclusionSize = 2000;
q = connectEM.generateQueriesFromAgglos(p, segmentMeta, agglosGtAxon2, outputFolder, options);
connectEM.debugNewQueries(segmentMeta, agglosGtAxon2, q, outputFolder);

% Do it on whole dataset for query number vs. agglomerateSize plots
options.writeTasksToFile = false;
options.queryBoundingBoxSize = 1000;
options.datasetBorderExclusionSize = 2000;
q = connectEM.generateQueriesFromAgglos(p, segmentMeta, agglos, outputFolder, options);
save([outputFolder 'allQueries.mat'], 'q');

