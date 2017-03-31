% This script was used for testing automated agglomeration & connectome generation
% Author: Manuel Berning <manuel.berning@brain.mpg.de>
% See +connectEM for additional functions used here

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
segmentMeta = addSegmentClassInformation(p, segmentMeta);
% Load synapse scores from SynEM
synScore = load([p.saveFolder 'globalSynScores.mat']);
synScore.isSynapse = connectEM.synScoresToSynEdges(graph, synScore);
toc;

display('Removing segments detected by heuristics & small & disconnected segments:');
tic;
[graphCut, excClasses, excNames, excSize] = connectEM.cutGraph(p, graph, segmentMeta, borderMeta, outputFolder, 100, 1000);
toc;

display('Writing skeletons and mappings of excluded components:');
tic;
for i=1:length(excClasses)
    connectEM.generateSkeletonFromAgglo(graph.edges, segmentMeta.point, ...
        excClasses{i}, ...
        strseq([excNames{i} '_'], 1:length(excClasses{i})), ...
        outputFolder, segmentMeta.maxSegId);
end
toc;

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
sizeThesholdDendrite = 1e6;
display('Performing agglomeration on dendrite subgraph:');
tic;
[dendrites, dendriteSize] = partitionSortAndKeepOnlyLarge(graphCutDendrites, segmentMeta, probThresholdDendrite, sizeThresholdDendrite);
toc;

probThresholdAxon = 0.90;
sizeThresholdAxon = 1e6;
display('Performing agglomeration on axon subgraph:');
tic;
[axons, axonsSize] = partitionSortAndKeepOnlyLarge(graphCutAxons, segmentMeta, probThresholdAxon, sizeThresholdAxon);
toc;

display('Removing components bordering on vessel or nuclei segments from axon graph:');
tic;
axonNeighbours = cellfun(@(x)unique(cat(2, graph.neighbours{x})), axons, 'uni', 0);
% Find fraction of neighbours for each axon component that are vessel or nuclei segments
vesselIdAll = cat(1, vessel{:});
vesselNeighbours = cellfun(@(x)sum(ismember(x, vesselIdAll)), axonNeighbours);
vesselFraction = vesselNeighbours ./ cellfun(@numel, axonNeighbours);
nucleiIdAll = cat(1, nuclei{:});
nucleiNeighbours = cellfun(@(x)sum(ismember(x, nucleiIdAll)), axonNeighbours);
nucleiFraction = nucleiNeighbours ./ cellfun(@numel, axonNeighbours);
clear axonNeighbours vesselIdAll vesselNeighbours nucleiIdAll nucleiNeighbours;
% Endothelial cells border on vessel, not axons
endo = axons(vesselFraction > 0);
% ER is part of soma, not axons
er = axons(nucleiFraction > 0 & vesselFraction == 0);
% Keep the rest
axons = axons(nucleiFraction == 0 & vesselFraction == 0); 
toc;

display('Attaching ER to dendrite class: ');
tic;
dendritesWithER = connectEM.attachER(graph, dendrites, er);
toc;

% Attach spines to dendrite class
display('Attaching spines to dendrite class: ');
tic;
[dendritesWithSpines, spinePaths] = connectEM.attachSpines(graph, segmentMeta, dendritesWithER, 10);
toc;

display('Writing skeletons for debugging the process:');
tic;
connectEM.generateSkeletonFromAgglo(dendriteEdges, segmentMeta.point, ...
    dendrites, ...
    strseq(['dendrites_' num2str(thresholdDendrite) '_'], 1:length(dendrites)), ...
    outputFolder, segmentMeta.maxSegId);
connectEM.generateSkeletonFromAgglo(graph.edges, segmentMeta.point, ...
    dendritesWithSpines, ...
    strseq(['dendritesWithAttachments_' num2str(thresholdDendrite) '_'], 1:length(dendritesWithSpines)), ...
    outputFolder, segmentMeta.maxSegId);
connectEM.generateSkeletonFromAgglo(axonEdges, segmentMeta.point, ...
    axons, ...
    strseq(['axons_' num2str(thresholdAxon) '_'], 1:length(axons)), ...
    outputFolder, segmentMeta.maxSegId);
% Probably not needed after refactoring as endo should already be excluded
connectEM.generateSkeletonFromAgglo(graph.edges, segmentMeta.point, ...
    endo, ...
    strseq('endo_', 1:length(endo)), ...
    outputFolder, segmentMeta.maxSegId);
% Here we need to check whether new probabilities still misclassify ER segments
connectEM.generateSkeletonFromAgglo(graph.edges, segmentMeta.point, ...
    er, ...
    strseq('er_', 1:length(er)), ...
    outputFolder, segmentMeta.maxSegId);
toc;

display('Display collected volume, save everything (except graph, which will be ignored due to size):');
tic;
% Take all classes (except removed small and disconnected segments as they do not have any good meaning)
eqClasses = cat(1, excClasses{1:6}, dendritesWithSpines, axons);
voxelCollected = sum(segmentMeta.voxelCount(cat(1, eqClasses{:})));
voxelTotal = sum(segmentMeta.voxelCount);
display(['Fraction of total (foreground) voxel collected: ' num2str(voxelCollected./voxelTotal, '%3.2f')]);
save([outputFolder 'agglo.mat']);
toc;

