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
% typeEM information ignorenum2str(p.raw.voxelSize(1))
%segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);

borderMeta = load([p.saveFolder 'globalBorder.mat'], 'borderSize', 'borderCoM');

synScore = load([p.saveFolder 'globalSynapseScores.mat']);
synScore.isSynapse = connectEM.synScoresToSynEdges(graph, synScore);
toc;

graphCut = connectEM.cutGraphSimple(p, graph, segmentMeta, borderMeta, 100, 1000);

% Lets stick with 99% for now as we have 'large enough' components
probThreshold = 0.99;
sizeThreshold = 1e6;
display(['Performing agglomeration on graph with thr prob:' num2str(probThreshold), ' size:' num2str(sizeThreshold)]);
tic;
[agglos, agglosSize, agglosEdges] = connectEM.partitionSortAndKeepOnlyLarge(graphCut, segmentMeta, probThreshold, sizeThreshold);
toc;

[agglosSizeSorted,idxSort] = sort(agglosSize,'descend');
agglosSorted= agglos(idxSort);

agglosOut = agglosSorted(1:100);
display('Writing skeletons for debugging the process:');

parameters.experiment.name= p.experimentName;
parameters.scale.x = num2str(p.raw.voxelSize(1));
parameters.scale.y = num2str(p.raw.voxelSize(2));
parameters.scale.z = num2str(p.raw.voxelSize(3));
parameters.offset.x = '0';
parameters.offset.y = '0';
parameters.offset.z = '0';

outputFolderSub = fullfile(outputFolder,['prob_' num2str(probThreshold) '_size' num2str(sizeThreshold)]);
mkdir(outputFolderSub)
tic;
Superagglos.skeletonFromAgglo(graph.edges, segmentMeta, ...
    agglosOut, 'agglos', outputFolderSub, parameters);
toc;

Util.save(fullfile(outputFolderSub,'agglos.mat'),agglos, agglosSize)

%{
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
%}
end
