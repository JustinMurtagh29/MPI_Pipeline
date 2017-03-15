function agglomerate(p)
% This script was used for testing automated agglomeration + focus flight queries for reconstruction of axons in 07x2
% Author: Manuel Berning <manuel.berning@brain.mpg.de>
% See agglomeration subfolder for additional functions used here

%% Start by loading parameter file

display('Loading data:');
tic;
%% Define where to store results
outputFolder = [p.saveFolder,'agglomeration/'];
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
%% Load graph and edge and segment (segments) based statistics
% Load global graph representation
graph = load([p.saveFolder 'graph.mat']);
% % Load information about edges
% borderMeta = load([p.saveFolder 'globalBorder.mat']);
% Load meta information of segments
segmentMeta = load([p.saveFolder 'segmentMeta.mat']);
%% Load synapse scores from SynEM
% synScore = load([p.saveFolder 'globalSynScores.mat']);
% synScore.isSynapse = connectEM.synScoresToSynEdges(graph, synScore);
toc;
% 
% display('Removing segments detected by heuristics:');
% tic;
% % Segments classified by heuristics
% load([p.saveFolder 'heuristicResult.mat']);
% excludedIds = segIds(vesselScore > 0.5 | myelinScore > 0.5 | nucleiScore > 0.5);
% % Remove only in edges
% keepIdx = ~any(ismember(graph.edges, excludedIds),2);
% graphCut.edges = graph.edges(keepIdx,:);
% graphCut.prob = graph.prob(keepIdx);
% toc;
graphCut = graph;

rng default;
threshold=0.99;%:-0.01:0.91;
initialPartition = cell(length(threshold),1);
sizePartition = initialPartition;
for t=1:length(threshold)
    display('Performing agglomeration:');
    tic;
    % Agglomerate segments using only GP probabilties
    [initialPartition{t}, remainingEdges] = connectEM.partitionWholeDataset(graphCut, threshold(t)); 
    sizePartition{t} = cellfun(@(x)sum(segmentMeta.voxelCount(x)), initialPartition{t});
    toc;
    % Probably make own function if needed again
    display('Writing skeletons for debugging the process:');
    tic;
    [~, idx] = sort(sizePartition{t}, 'descend');
    parameters.experiment.name = p.experimentName;
    parameters.scale = struct('x',p.raw.voxelSize(1),'y',p.raw.voxelSize(2),'z',p.raw.voxelSize(3));
    connectEM.generateSkeletonFromAgglo(remainingEdges, segmentMeta.point', ...
        initialPartition{t}(idx(1:1000)), ...
        strseq(['largestComponents' num2str(threshold(t)) '_'], 1:1000), ...
        outputFolder, segmentMeta.maxSegId,parameters);
    %idx = randi(numel(initialPartition{t}),100,1);
    %connectEM.generateSkeletonFromAgglo(remainingEdges, segmentMeta.point', ...
    %    initialPartition{t}(idx), ...
    %    strseq(['randomComponents' num2str(threshold(t)) '_'], 1:100), ...
    %    outputFolder, segmentMeta.maxSegId);
    toc;
end
voxelCount = segmentMeta.voxelCount;
maxSegId = segmentMeta.maxSegId;
save([outputFolder 'initialAgglo.mat'], 'initialPartition', 'sizePartition', 'voxelCount', 'maxSegId');

