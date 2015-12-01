% Load needed data files
load('/gaba/u/mberning/results/pipeline/20151111T183414/allParameter.mat');
load([p.saveFolder 'graph.mat']);
load([p.saveFolder 'CoM.mat']);
seedCube = p.local(6,9,9).bboxSmall;

% Restrict to some region
[graphR, comR, segIdsR, segIdsB] = restrictSGtoRegion(p, graph, com, seedCube);
% Cluster or threshold or similar?

% Find which edges to annotate

% Write to knowledge database


