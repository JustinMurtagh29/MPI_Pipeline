% This script was used for testing automated agglomeration & connectome generation
% Author: Manuel Berning <manuel.berning@brain.mpg.de>
% See +connectEM for additional functions used here

clearvars -except graph borderMeta segmentMeta synScore graphCut excClasses excNames excSize;
%% Start by loading parameter file
% Load parameter from newest pipeline run
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
% To keep workspace clean here (and we are gonna have a bunch of stuff anyway) remove parameter for training (pT)
clear pT;
%% Define where to store results
outputFolder = ['/gaba/scratch/mberning/' datestr(clock,30) '_agglomeration/'];
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%{
display('Loading data:');
tic;
% Load global graph representation
graph = load([p.saveFolder 'graph.mat'], 'prob', 'edges', 'neighbours', 'neighProb');
% Load information about edges
borderMeta = load([p.saveFolder 'globalBorder.mat'], 'borderSize', 'borderCoM');
% Load meta information of segments
segmentMeta = load([p.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId');
segmentMeta.point = segmentMeta.point';
% Load and preprocess segment class predictions on single segments from Alessandro
segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);
% Load synapse scores from SynEM
synScore = load([p.saveFolder 'globalSynScores.mat']);
synScore.isSynapse = connectEM.synScoresToSynEdges(graph, synScore);
toc;

display('Removing segments detected by heuristics & small & disconnected segments:');
tic;
[graphCut, excClasses, excNames, excSize] = connectEM.cutGraph(p, graph, segmentMeta, borderMeta, 100, 1000);
toc;

display('Writing skeletons (probably not useable due to not being CC) and mappings of excluded components:');
tic;
connectEM.skeletonFromAgglo(graph.edges, segmentMeta, ...
    excClasses(1:6), 'heuristics', outputFolder);
connectEM.skeletonFromAgglo(graph.edges, segmentMeta, ...
    excClasses(7:8), 'smallAndDisconnected', outputFolder);
toc;
%}

display('Generating subgraphs for axon and dendrite agglomeration:');
tic;
% Dendrites first
idx = all(ismember(graphCut.edges, find(segmentMeta.isDendrite)), 2);
graphCutDendrites.edges = graphCut.edges(idx,:);
graphCutDendrites.prob = graphCut.prob(idx);
% Then axon
idx = all(ismember(graphCut.edges, find(segmentMeta.isAxon)), 2);
graphCutAxons.edges = graphCut.edges(idx,:);
graphCutAxons.prob = graphCut.prob(idx);
clear idx;
toc;

% Lets stick with 99% for now as we have 'large enough' components
probThresholdDendrite = 0.99;
sizeThresholdDendrite = 1e6;
display('Performing agglomeration on dendrite subgraph:');
tic;
[dendrites, dendriteSize, dendriteEdges] = connectEM.partitionSortAndKeepOnlyLarge(graphCutDendrites, segmentMeta, probThresholdDendrite, sizeThresholdDendrite);
toc;

probThresholdAxon = 0.90;
sizeThresholdAxon = 1e6;
display('Performing agglomeration on axon subgraph:');
tic;
[axons, axonsSize, axonEdges] = connectEM.partitionSortAndKeepOnlyLarge(graphCutAxons, segmentMeta, probThresholdAxon, sizeThresholdAxon);
toc;

display('Reassigning ER from axon to dendrite class: ');
tic;
[dendritesAfterEr, axonsAfterEr, er] = connectEM.extractAndTransferER(graph, dendrites, axons);
toc;

display('Garbage collection');
tic;
[axonsFinal, dendritesFinal] = connectEM.garbageCollection(graph, segmentMeta, axonsAfterEr, dendritesAfterEr, excClasses(1:6));
toc;

display('Attaching spines to dendrite class: ');
tic;
[dendritesFinalWithSpines, spinePaths, comment] = connectEM.attachSpines(graph, segmentMeta, dendritesFinal, 10);
toc;

display('Writing skeletons for debugging the process:');
tic;
connectEM.skeletonFromAgglo(graph.edges, segmentMeta, ...
    er, 'er', outputFolder);
rng default;
randIdx = randi(numel(spinePaths), 100);
thisSpinePath = spinePaths(randIdx);
thisTreeName = cellfun(@(x,y)[num2str(y) '_' x], comment(randIdx), num2cell(1:100));
thisNodes = cellfun(@(x)segmentMeta.point(x,:), thisSpinePath);
connectEM.generateSkeletonFromNodes([outputFolder 'spineTestSet.nml'], thisNodes, thisTreeName, {}, true);
connectEM.skeletonFromAgglo(graph.edges, segmentMeta, ...
    er, 'er', outputFolder);
connectEM.skeletonFromAgglo(graph.edges, segmentMeta, ...
    dendritesFinal, 'dendrites', outputFolder);
connectEM.skeletonFromAgglo(graph.edges, segmentMeta, ...
    dendritesFinal, 'dendrites', outputFolder);
connectEM.skeletonFromAgglo(graph.edges, segmentMeta, ...
    dendritesFinalWithSpines, 'dendritesWithSpines', outputFolder);
connectEM.skeletonFromAgglo(axonEdges, segmentMeta, ...
    axonsFinal, 'axons', outputFolder);
toc;

display('Display collected volume, save everything (except graph, which will be ignored due to size):');
tic;
% Take all classes (except removed small and disconnected segments as they do not have any good meaning)
collectedSegments = false(segmentMeta.maxSegId, 1);
collectedSegments(cat(1, excClasses{1:6}, dendritesWithSpines{:}, axons{:})) = true;
voxelCollected = sum(segmentMeta.voxelCount(collectedSegments & segmentMeta.isDendrite));
voxelTotal = sum(segmentMeta.voxelCount(segmentMeta.isDendrite));
display(['Fraction of total dendrite voxel collected: ' num2str(voxelCollected./voxelTotal, '%3.2f')]);
voxelCollected = sum(segmentMeta.voxelCount(collectedSegments & segmentMeta.isAxon));
voxelTotal = sum(segmentMeta.voxelCount(segmentMeta.isAxon));
display(['Fraction of total axon voxel collected: ' num2str(voxelCollected./voxelTotal, '%3.2f')]);
save([outputFolder 'agglo.mat']);
toc;

display('Generating isosurfaces of final axon and dendrite components:');
tic;

toc;

