
% load data
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat', 'voxelCount', 'centroid', 'box', 'maxSegId', 'cubeIdx', 'point');
a = load('/gaba/scratch/kboerg/cycleAgglos_forQueries_forMB.mat');
b = load('/gaba/scratch/kboerg/directCycle/cycle003_for_mb.mat');
agglos = b.axonsFinal(a.y_zz.axon1.foundAgglomerates_col{2});
clear a b;
% calculate queries and write to file on wk
outputFolder = '/gaba/scratch/mberning/flightQueries/';
q = connectEM.generateQueriesFromAgglos(p, segmentMeta, agglos, outputFolder);
connectEM.debugNewQueries(segmentMeta, agglos, q, outputFolder);

