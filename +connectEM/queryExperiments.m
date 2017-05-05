
% load data
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat', 'voxelCount', 'centroid', 'box', 'maxSegId', 'cubeIdx', 'point');
a = load('/gaba/scratch/kboerg/cycleAgglos_forQueries_forMB.mat');
agglos = a.y_zz.axon1.foundAgglomerates_col;
% calculate queries and write to file on wk
outputFolder = '/gaba/scratch/mberning/flightQueries/';
q = connectEM.generateQueriesFromAgglos(p, segmentMeta, agglos, outputFolder);
connectEM.debugAgglomerationWithNewQueries(segmentMeta, agglos, q, outputFolder);
