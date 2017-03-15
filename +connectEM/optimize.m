% Load parameter struct, graph & segment information needed
load /gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat;
graph = load([p.saveFolder 'graph.mat'], 'edges', 'prob');
segMeta = load([p.saveFolder 'segmentMeta.mat'], 'point', 'voxelCount', 'maxSegId');
segMeta.point = segMeta.point';

% Keep only segments larger than 100 voxel
segmentsToKeep = find(segMeta.voxelCount > 100);
idx = all(ismember(graph.edges, segmentsToKeep),2);
% Also keep continuos numbering, should make steps below easier for now
graphCut.edges = connectEM.changem(graph.edges(idx,:), 1:numel(segmentsToKeep), segmentsToKeep);
graphCut.prob = graph.prob(idx);
maxSegId = numel(segmentsToKeep);
clear idx;

% Generate better representation
probM = accumarray(double(graphCut.edges), graphCut.prob, ...
    double(repmat(maxSegId, 1, 2)), @max, 0, true);
% Inital guess
x0 = probM > 0.99;

% Evaluate objective function at starting conditions
tic; value = connectEM.aggloObjective(probM, x0); toc;

% Try some global optimization approaches
opts = optimoptions(@fmincon,'Algorithm','interior-point');
problem = createOptimProblem('fmincon', whatelse);
