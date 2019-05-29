function agglomerate(p);
% This script was used for testing automated agglomeration & connectome generation
% Author: 
%           Manuel Berning <manuel.berning@brain.mpg.de>
% Modified by:
%           Sahil Loomba <sahil.loomba@brain.mpg.de>
% See +connectEM for additional functions used here

outputFolder = [p.saveFolder datestr(clock,30) '_agglomeration/'];                                                                                                           
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

display('loading dependent variables...');
tic;
m = load([p.saveFolder 'allParameterWithSynapses.mat']);
p = m.p;

graph = load([p.saveFolder 'graph.mat']);

segmentMeta = load([p.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId');
segmentMeta.point = segmentMeta.point';
segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);

borderMeta = load([p.saveFolder 'globalBorder.mat'], 'borderSize', 'borderCoM');

synScore = load([p.saveFolder 'globalSynapseScores.mat'])
synScore.isSynapse = connectEM.synScoresToSynEdges(graph, synScore);
toc;

forceKeepEdges = []; 
% added by SL temp fix
segmentMeta.cubeIdx = [];
m = load([p.saveFolder 'heuristicResult.mat']);
heuristics.mapping = m.segIds;

graphCut = connectEM.cutGraph(p, graph, segmentMeta, borderMeta, heuristics, 100, 1000,forceKeepEdges);

%{
display('Removing segments detected by heuristics & small & disconnected segments:');
tic;
[graphCut, excClasses, excNames, excSize] = connectEM.cutGraph(p, graph, segmentMeta, borderMeta, 100, 1000);
toc;

display('Writing skeletons (probably not useable due to not being CC) and mappings of excluded components:');
tic;
Superagglos.skeletonFromAgglo(graph.edges, segmentMeta, ...
    excClasses(1:6), 'heuristics', outputFolder);
Superagglos.skeletonFromAgglo(graph.edges, segmentMeta, ...
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
[dendritesFinalWithSpines, spinePaths, comment] = connectEM.attachSpines(graph, segmentMeta, dendritesFinal, axonsFinal, 10);
toc;

display('Writing skeletons for debugging the process:');
tic;
Superagglos.skeletonFromAgglo(graph.edges, segmentMeta, ...
    er, 'er', outputFolder);
Superagglos.skeletonFromAgglo(graph.edges, segmentMeta, ...
    er, 'er', outputFolder);
Superagglos.skeletonFromAgglo(graph.edges, segmentMeta, ...
    dendritesFinal, 'dendrites', outputFolder);
Superagglos.skeletonFromAgglo(graph.edges, segmentMeta, ...
    dendritesFinal, 'dendrites', outputFolder);
Superagglos.skeletonFromAgglo(graph.edges, segmentMeta, ...
    dendritesFinalWithSpines, 'dendritesWithSpines', outputFolder);
Superagglos.skeletonFromAgglo(axonEdges, segmentMeta, ...
    axonsFinal, 'axons', outputFolder);
toc;

display('Display collected volume, save everything (except graph, which will be ignored due to size):');
tic;
% Take all classes (except removed small and disconnected segments as they do not have any good meaning)
collectedSegments = false(segmentMeta.maxSegId, 1);
collectedSegments(cat(1, excClasses{1:6}, dendritesFinalWithSpines{:}, axons{:})) = true;
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
job1 = connectEM.buildIsosurfaceOfAggloStart(p, outputFolder, dendritesFinalWithSpines, 'dendritesFinalWithSpines');
job2 = connectEM.buildIsosurfaceOfAggloStart(p, outputFolder, axonsFinal, 'axonsFinal');
job3 = connectEM.buildIsosurfaceOfAggloStart(p, outputFolder, excClasses(1:6)', 'heuristics');
toc;

display('Write 100 random dendritic components and 100 spine paths for error annotation:');
tic;
rng default;
randIdx = randperm(length(dendritesFinalWithSpines), 100);
Superagglos.skeletonFromAgglo(graph.edges, segmentMeta, ...
    dendritesFinal(randIdx), 'dendritesForErrorAnnotation', outputFolder);
rng default;
randIdx = randperm(numel(spinePaths), 100);
thisSpinePath = spinePaths(randIdx);
thisTreeName = cellfun(@(x,y)[num2str(y) '_' x], comment(randIdx), num2cell(1:100)', 'uni', 0);
thisNodes = cellfun(@(x)segmentMeta.point(x,:), thisSpinePath, 'uni', 0);
connectEM.generateSkeletonFromNodes([outputFolder 'spinesForErrorAnnotation.nml'], thisNodes, thisTreeName, {}, true);
toc;

% Split up nuclei and render individually
[nuclei, nucleiSize] = connectEM.findCCaccordingToGraph(graph, excClasses{3}, segmentMeta);
% Determine which cell each nucleus belongs to
resortedNuclei = connectEM.attachNuclei(graph, nuclei, dendrites);
job4 = connectEM.buildIsosurfaceOfAggloStart(p, outputFolder, resortedNuclei', 'nuclei');
myelin = connectEM.findCCaccordingToGraph(graph, excClasses{5}, segmentMeta, 1e5);
job5 = connectEM.buildIsosurfaceOfAggloStart(p, outputFolder, myelin', 'myelin');

% New axon agglomeration, old one was shit!
% Then axon
tic;
idx = all(ismember(graphCut.edges, find(segmentMeta.axonProb > .9)), 2);
graphCutAxons.edges = graphCut.edges(idx,:);
graphCutAxons.prob = graphCut.prob(idx);
probThresholdAxon = 0.9;
sizeThresholdAxon = 5e5;
[axonsNew, axonsSize, axonEdges] = connectEM.partitionSortAndKeepOnlyLarge(graphCutAxons, segmentMeta, probThresholdAxon, sizeThresholdAxon);
Superagglos.skeletonFromAgglo(axonEdges, segmentMeta, ...
        axonsNew, 'axonsNew', outputFolder);
toc;

% Select "good" axons
[~, ~, latent] = cellfun(@(x)princomp(segmentMeta.point(x,:)), axonsNew);
explainedVariance = cellfun(@(x)cumsum(x)./sum(x), latent, 'uni', 0);
explainedVariance1 = cellfun(@(x)x(1), explainedVariance, 'uni', 0);
axonsGood = axonsNew(explainedVariance1 > 0.95);
Superagglos.skeletonFromAgglo(axonEdges, segmentMeta, ...
        axonsGood, 'axonsGood', outputFolder);
jobAxons = connectEM.buildIsosurfaceOfAggloStart(p, outputFolder, axonsGood', 'axonsGood');

% Write spines in 5 micron^3
spineSegments = find(segmentMeta.isSpine);
spinePosition = segmentMeta.point(spineSegments,:);
bbox_wk = [2800, 4267, 1712, 445, 445, 179];
bbox = Util.convertWebknossosToMatlabBbox(bbox_wk);
idx = all(bsxfun(@minus, spinePosition, bbox(:,1)') > 0,2) & all(bsxfun(@minus, spinePosition, bbox(:,2)') < 0,2);
connectEM.generateSkeletonFromNodes([outputFolder 'spinesBbox.nml'], spinePosition(idx), strseq('spines', 1:sum(idx), {}, true));
end
