%function agglomerate(p);
% This script was used for testing automated agglomeration & connectome generation
% Author: 
%           Manuel Berning <manuel.berning@brain.mpg.de>
% Modified by:
%           Sahil Loomba <sahil.loomba@brain.mpg.de>
% See +connectEM for additional functions used here

info = Util.runInfo();

outputFolder = [p.saveFolder datestr(clock,30) '_agglomeration/'];
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

display('loading dependent variables...');
m = load([p.saveFolder 'allParameterWithSynapses.mat']);
p = m.p;

graph = load([p.saveFolder 'graph.mat']);

segmentMeta = load([p.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId');
segmentMeta.point = segmentMeta.point';
% typeEM information ignorenum2str(p.raw.voxelSize(1))
segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta); % SL: spine not added yet
borderMeta = load([p.saveFolder 'globalBorder.mat'], 'borderSize', 'borderCoM');
%synScore = load([p.saveFolder 'globalSynapseScores.mat']);
%synScore.isSynapse = connectEM.synScoresToSynEdges(graph, synScore);

maxSegId = Seg.Global.getMaxSegId(p);
borderSizeThr = 50;
segmentSizeThr = 100;
display(['Cut graph at border size:' num2str(borderSizeThr), 'segment size:' num2str(segmentSizeThr)]);
[graphCut, forceKeepEdges] = connectEM.cutGraphSimple(p, graph, segmentMeta, borderMeta, borderSizeThr, segmentSizeThr);

Util.log('Update segment type prob for corr edges to the max of the edge:')
corrIdx = find(isnan(graph.borderIdx));
for i=1:size(corrIdx,1)
    curEdges = graph.edges(corrIdx(i),:);
    segmentMeta.dendriteProb(curEdges) = max(segmentMeta.dendriteProb(curEdges));
    segmentMeta.axonProb(curEdges) = max(segmentMeta.axonProb(curEdges));
    segmentMeta.gliaProb(curEdges) = max(segmentMeta.gliaProb(curEdges));
end

Util.log('Setting type prob thresholds for classes:')
segmentMeta.isDendrite = segmentMeta.dendriteProb > 0.85 & segmentMeta.axonProb < 0.50;
segmentMeta.isAxon = segmentMeta.axonProb >= 0.90;

Util.log('Generating subgraphs for axon and dendrite agglomeration:');
% Dendrites first: both edge seg belong to dend
tempClass = ismember(graphCut.edges, find(segmentMeta.isDendrite));
idxDend = all(tempClass, 2);
forceKeepEdgesDend = false(size(forceKeepEdges,1),1);
forceKeepEdgesDend(forceKeepEdges) = tempClass(forceKeepEdges,1) | tempClass(forceKeepEdges,2); % keep those corr edges which has either seg as dendrite
idx = idxDend | forceKeepEdgesDend;
graphCutDendrites.edges = graphCut.edges(idx,:);
graphCutDendrites.prob = graphCut.prob(idx);

% Then axon
tempClass = ismember(graphCut.edges, find(segmentMeta.isAxon));
idxAxon = all(tempClass, 2);
forceKeepEdgesAxon = false(size(forceKeepEdges,1),1);
forceKeepEdgesAxon(forceKeepEdges) = tempClass(forceKeepEdges,1) | tempClass(forceKeepEdges,2); % keep those corr edges which has either seg as axon
idx = idxAxon | forceKeepEdgesAxon;
graphCutAxons.edges = graphCut.edges(idx,:);
graphCutAxons.prob = graphCut.prob(idx);
clear idx;

%% Lets stick with 99% for now as we have 'large enough' components
probThresholdDend = 0.99;
sizeThresholdDend = 1e3;
Util.log(['Performing agglomeration on dendrite subgraph with thr prob:' num2str(probThresholdDend), 'agglo size:' num2str(sizeThresholdDend)]);
[dendrites, dendriteSize, dendriteEdges] = connectEM.partitionSortAndKeepOnlyLarge(graphCutDendrites, segmentMeta,...
                                                         probThresholdDend, sizeThresholdDend);
[dendriteSizeSorted,idxSort] = sort(dendriteSize,'descend');
dendritesSorted = dendrites(idxSort);
Util.log('Calculating dendrite path lengths...')
dendritePathLengths = cellfun(@(x) connectEM.pathLengthOfAgglo(segmentMeta, x, [11.24, 11.24, 30]), dendritesSorted,'uni',0); % um
dendritePathLengths = cell2mat(dendritePathLengths);
outputFolderSub = fullfile(outputFolder,['dendrites_border_' num2str(borderSizeThr) ...
                            'seg_' num2str(segmentSizeThr) ...
                            'prob_' num2str(probThresholdDend) ...
                            'size_' num2str(sizeThresholdDend)]);
mkdir(outputFolderSub);
Util.save(fullfile(outputFolderSub,'dendrites.mat'),dendritesSorted, dendriteSizeSorted, dendritePathLengths, info)

Util.log('Write new segmentation based on agglos')
segOut = struct;
segOut.root = fullfile(p.saveFolder, 'aggloState', ...
    [datestr(clock,30) '_segForDendAgglos'], '1');
mkdir(segOut.root)
dataType = 'uint32';
wkwInit('new',segOut.root,32, 32, dataType, 1);
segOut.backend = 'wkwrap';
%segOut.prefix = p.seg.prefix;
agglosSorted = dendritesSorted;
mapping = connectEM.createLookup(segmentMeta, agglosSorted);
Seg.Global.applyMappingToSegmentation(p, mapping, segOut);
thisBBox = [1, 1, 1; (ceil(p.bbox(:, 2) ./ 1024) .* 1024)']';
createResolutionPyramid(segOut, thisBBox, [], true);

%{
% Export skeletons
agglosOut = dendritesSorted(1:100);
display('Writing skeletons for debugging the process:');
parameters.experiment.name= p.experimentName;
parameters.scale.x = num2str(p.raw.voxelSize(1));
parameters.scale.y = num2str(p.raw.voxelSize(2));
parameters.scale.z = num2str(p.raw.voxelSize(3));
parameters.offset.x = '0';
parameters.offset.y = '0';
parameters.offset.z = '0';
Superagglos.skeletonFromAgglo(graphCutDendrites.edges, segmentMeta, ...
    agglosOut, 'dendrites', outputFolderSub, parameters);
%}

%% repeat for axons
probThresholdAxon = 0.99;
sizeThresholdAxon = 1e3;
Util.log(['Performing agglomeration on axon subgraph with thr prob:' num2str(probThresholdAxon), 'agglo size:' num2str(sizeThresholdAxon)]);
[axons, axonSize, axonEdges] = connectEM.partitionSortAndKeepOnlyLarge(graphCutAxons, segmentMeta,...
                                                         probThresholdAxon, sizeThresholdAxon);
[axonSizeSorted,idxSort] = sort(axonSize,'descend');
axonsSorted = axons(idxSort);
Util.log('Calculating axon path lengths...')
axonPathLengths = cellfun(@(x) connectEM.pathLengthOfAgglo(segmentMeta, x, [11.24, 11.24, 30]), axonsSorted,'uni',0);% um
axonPathLengths = cell2mat(axonPathLengths);
outputFolderSub = fullfile(outputFolder,['axons_border_' num2str(borderSizeThr) ...
                            'seg_' num2str(segmentSizeThr) ...
                            'prob_' num2str(probThresholdAxon) ...
                            'size_' num2str(sizeThresholdAxon)]);
mkdir(outputFolderSub);
Util.save(fullfile(outputFolderSub,'axons.mat'),axonsSorted, axonSizeSorted, axonPathLengths, info)

Util.log('Write new segmentation based on agglos')
segOut = struct;
segOut.root = fullfile(p.saveFolder, 'aggloState', ...
    [datestr(clock,30) '_segForAxonAgglos'], '1');
mkdir(segOut.root)
dataType = 'uint32';
wkwInit('new',segOut.root,32, 32, dataType, 1);
segOut.backend = 'wkwrap';
%segOut.prefix = p.seg.prefix;
agglosSorted = axonsSorted;
mapping = connectEM.createLookup(segmentMeta, agglosSorted);
Seg.Global.applyMappingToSegmentation(p, mapping, segOut);
thisBBox = [1, 1, 1; (ceil(p.bbox(:, 2) ./ 1024) .* 1024)']';
createResolutionPyramid(segOut, thisBBox, [], true);


%{
% Export skeletons
agglosOut = axonsSorted(1:100);
display('Writing skeletons for debugging the process:');
parameters.experiment.name= p.experimentName;
parameters.scale.x = num2str(p.raw.voxelSize(1));
parameters.scale.y = num2str(p.raw.voxelSize(2));
parameters.scale.z = num2str(p.raw.voxelSize(3));
parameters.offset.x = '0';
parameters.offset.y = '0';
parameters.offset.z = '0';
Superagglos.skeletonFromAgglo(graphCutAxons.edges, segmentMeta, ...
%}
   
%{
%% plot agglo length statistics
Util.log('Now plotting path lengths...')
fig = figure;
fig.Color = 'white';
subplot(1,2,1)
[n, xout] = hist(dendritePathLengths,100);
bar(xout, n, 'barwidth', 1, 'basevalue', 1);
set(gca,'YScale','log')
ylim([0, 1e7]);
Util.setPlotDefault(gca,0,0);
title('dendrites')
xlabel('Path length (\mum)');
ylabel('Frequency (log)')
subplot(1,2,2)
[n, xout] = hist(axonPathLengths,100);
bar(xout, n, 'barwidth', 1, 'basevalue', 1);
set(gca,'YScale','log')
ylim([0, 1e7]);
Util.setPlotDefault(gca,0,0);
title('axons')
xlabel('Path length (\mum)');
ylabel('Frequency (log)')
saveas(gcf,fullfile(outputFolder,'agglo_path_lengths.png'))
close all
%}

%{
display('Garbage collection');
tic;
[axonsFinal, dendritesFinal] = connectEM.garbageCollection(graph, segmentMeta, axonsAfterEr, dendritesAfterEr, excClasses(1:6));
toc;
%}

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
%end
