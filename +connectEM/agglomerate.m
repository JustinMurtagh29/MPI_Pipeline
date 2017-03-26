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
%% Load synapse scores from SynEM
synScore = load([p.saveFolder 'globalSynScores.mat']);
synScore.isSynapse = connectEM.synScoresToSynEdges(graph, synScore);
toc;

display('Removing segments detected by heuristics:');
tic;
[graphCut, eqClassesExcluded] = connectEM.cutGraph(p, segmentMeta, borderMeta, outputFolder, 100, 1000); 
[vessel, sizeVessel] = connectEM.findCCaccordingToGraph(graph, eqClassesExcluded{1}, segmentMeta);
[nuclei, sizeNuclei] = connectEM.findCCaccordingToGraph(graph, eqClassesExcluded{2}, segmentMeta);
[myelin, sizeMyelin] = connectEM.findCCaccordingToGraph(graph, eqClassesExcluded{3}, segmentMeta);
toc;

% Lets stick with 99% for now as we have 'large enough' components
rng default;
threshold = 0.99;
display('Performing agglomeration:');
tic;
% Agglomerate segments using only GP probabilties
[initialPartition, remainingEdges] = connectEM.partitionWholeDataset(graphCut, threshold); 
sizePartition = cellfun(@(x)sum(segmentMeta.voxelCount(x)), initialPartition);
[sizePartition, idx] = sort(sizePartition, 'descend');
initialPartition = initialPartition(idx);
% Keep only agglomerates that have at least 1 million voxels
idx = sizePartition > 1e6; 
initialPartition = initialPartition(idx);
sizePartition = sizePartition(idx);
clear idx;
toc;

% Probably make own function if needed again
display('Writing skeletons for debugging the process:');
tic;
connectEM.generateSkeletonFromAgglo(remainingEdges, segmentMeta.point, ...
    initialPartition, ...
    strseq(['allComponents' num2str(threshold) '_'], 1:length(initialPartition)), ...
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
eqClasses = cat(1, initialPartition, vessel, nuclei, myelin);
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

