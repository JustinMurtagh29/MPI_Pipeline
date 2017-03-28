% This script was used for testing automated agglomeration + focus flight queries for reconstruction of axons in 07x2
% Author: Manuel Berning <manuel.berning@brain.mpg.de>
% See agglomeration subfolder for additional functions used here

%% Start by loading parameter file
% Load parameter from newest pipeline run
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
% To keep workspace clean here (and we are gonna have a bunch of stuff anyway) remove parameter for training (pT)
clear pT;

display('Loading data:');
tic;
%% Define where to store results
outputFolder = ['/gaba/scratch/mberning/' datestr(clock,30) '_agglomeration/'];
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
%% Load graph and edge and segment (segments) based statistics
% Load global graph representation
graph = load([p.saveFolder 'graph.mat'], 'prob', 'edges');
% Load information about edges
borderMeta = load([p.saveFolder 'globalBorder.mat'], 'borderSize', 'borderCoM');
% Load meta information of segments
segmentMeta = load([p.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId');
segmentMeta.point = segmentMeta.point';
% Load and preprocess segment class predictions on single segments from Alessandro
temp = load([p.saveFolder 'segmentPredictions.mat']);
% ... axon
segmentMeta.axonProb = zeros(segmentMeta.maxSegId, 1);
idx = ~isnan(temp.probs(:,2));
segmentMeta.axonProb(temp.segId(idx)) = temp.probs(idx,2);
segmentMeta.isAxon = segmentMeta.axonProb > 0.42;
% ... dendrite
segmentMeta.dendriteProb = zeros(segmentMeta.maxSegId, 1);
idx = ~isnan(temp.probs(:,3));
segmentMeta.dendriteProb(temp.segId(idx)) = temp.probs(idx,3);
segmentMeta.isDendrite = segmentMeta.dendriteProb > 0.20;
clear temp idx;
%% Load synapse scores from SynEM
synScore = load([p.saveFolder 'globalSynScores.mat']);
synScore.isSynapse = connectEM.synScoresToSynEdges(graph, synScore);
toc;

display('Removing segments detected by heuristics & small & disconnected segments:');
tic;
[graphCut, eqClassesExcluded] = connectEM.cutGraph(p, segmentMeta, borderMeta, outputFolder, 100, 1000);
[vessel, sizeVessel] = connectEM.findCCaccordingToGraph(graph, eqClassesExcluded{1}, segmentMeta);
[nuclei, sizeNuclei] = connectEM.findCCaccordingToGraph(graph, eqClassesExcluded{2}, segmentMeta);
[myelin, sizeMyelin] = connectEM.findCCaccordingToGraph(graph, eqClassesExcluded{3}, segmentMeta);
toc;

display('Generating subgraphs for axon and dendrite agglomeration:');
tic;
idx = all(ismember(graphCut.edges, find(segmentMeta.isDendrite)), 2);
graphCutDendrites.edges = graphCut.edges(idx,:);
graphCutDendrites.prob = graphCut.prob(idx);
idx = all(ismember(graphCut.edges, find(segmentMeta.isAxon)), 2);
graphCutAxons.edges = graphCut.edges(idx,:);
graphCutAxons.prob = graphCut.prob(idx);
clear idx;
toc;

% Lets stick with 99% for now as we have 'large enough' components
rng default;
thresholdDendrite = 0.99;
display('Performing agglomeration on dendrite subgraph:');
tic;
% Agglomerate segments using only GP probabilties
[dendrites, dendriteEdges] = connectEM.partitionWholeDataset(graphCutDendrites, thresholdDendrite); 
dendriteSize = cellfun(@(x)sum(segmentMeta.voxelCount(x)), dendrites);
[dendriteSize, idx] = sort(dendriteSize, 'descend');
dendrites = dendrites(idx);
% Keep only agglomerates that have at least 1 million voxels
idx = dendriteSize > 1e6; 
dendrites = dendrites(idx);
dendriteSize = dendriteSize(idx);
clear idx;
toc;

rng default;
thresholdAxon = 0.99;
display('Performing agglomeration on axon subgraph:');
tic;
% Agglomerate segments using only GP probabilties
[axons, axonEdges] = connectEM.partitionWholeDataset(graphCutAxons, thresholdAxon); 
axonSize = cellfun(@(x)sum(segmentMeta.voxelCount(x)), axons);
[axonSize, idx] = sort(axonSize, 'descend');
axons = axons(idx);
% Keep only agglomerates that have at least 1 million voxels
idx = axonSize > 1e6; 
axons = axons(idx);
axonSize = axonSize(idx);
clear idx;
toc;

display('Writing skeletons for debugging the process:');
tic;
connectEM.generateSkeletonFromAgglo(dendriteEdges, segmentMeta.point, ...
    dendrites, ...
    strseq(['dendrites_' num2str(thresholdDendrite) '_'], 1:length(dendrites)), ...
    outputFolder, segmentMeta.maxSegId);
connectEM.generateSkeletonFromAgglo(axonEdges, segmentMeta.point, ...
    axons, ...
    strseq(['axons_' num2str(thresholdAxon) '_'], 1:length(axons)), ...
    outputFolder, segmentMeta.maxSegId);
connectEM.generateSkeletonFromAgglo(graph.edges, segmentMeta.point, ...
    vessel, ...
    strseq('vessel_', 1:length(vessel)), ...
    outputFolder, segmentMeta.maxSegId);
connectEM.generateSkeletonFromAgglo(graph.edges, segmentMeta.point, ...
    nuclei, ...
    strseq('nuclei_', 1:length(nuclei)), ...
    outputFolder, segmentMeta.maxSegId);
connectEM.generateSkeletonFromAgglo(graph.edges, segmentMeta.point, ...
    myelin, ...
    strseq('myelin_', 1:length(myelin)), ...
    outputFolder, segmentMeta.maxSegId);
toc;

display('Display collected volume, write mapping of current state, save everything:');
tic;
eqClasses = cat(1, dendrites, axons, vessel, nuclei, myelin);
voxelCollected = sum(segmentMeta.voxelCount(cat(1, eqClasses{:})));
voxelTotal = sum(segmentMeta.voxelCount);
display(['Fraction of total (foreground) voxel collected: ' num2str(voxelCollected./voxelTotal, '%3.2f')]);
script = WK.makeMappingScript(segmentMeta.maxSegId, eqClasses);
fileHandle = fopen([outputFolder 'stateAfterInitialAgglo.txt'], 'w');
fwrite(fileHandle, script);
fclose(fileHandle);
clear fileHandle script;
save([outputFolder 'initialAgglo.mat']);
toc;

