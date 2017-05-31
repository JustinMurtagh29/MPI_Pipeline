% Load needed data
graph = load('/gaba/u/mberning/results/pipeline/20170217_ROI/graphNew.mat', 'edges', 'prob', 'borderIdx');
[graph.neighbours, neighboursIdx] = Graph.edges2Neighbors(graph.edges);
graph.neighProb = cellfun(@(x)graph.prob(x), neighboursIdx, 'uni', 0);
graph.neighBorderIdx = cellfun(@(x)graph.borderIdx(x), neighboursIdx, 'uni', 0);
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat');
segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);
borderMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalBorder.mat', 'borderSize', 'borderCoM');
globalSegmentPCA = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalSegmentPCA.mat', 'covMat');

% Options for algorithm
options.axonScore = 0.5;
options.bboxDist = 1000;
options.minSize = 100;
options.voxelSize = [11.24 11.24 28];

% Exclude all segments in cubes with catastrophic merger
[er, cm] = connectEM.getERcomponents();
excludedCubeIdx = unique(cellfun(@(x)mode(segmentMeta.cubeIdx(x)), cm));
excludedSegmentIdx = ismember(segmentMeta.cubeIdx, excludedCubeIdx) | ismember(1:segmentMeta.maxSegId, cat(1, er{:}))';
% Add single segments if not excluded due to catastrophic merger or ER 
startAgglo = load('/gaba/scratch/mberning/aggloGridSearch/search05_00564.mat', 'axonsFinal');
axons = cat(1, startAgglo.axonsFinal, ...
    num2cell(setdiff(find(segmentMeta.axonProb > options.axonScore & ~excludedSegmentIdx), cell2mat(startAgglo.axonsFinal))));
clear startAgglo excluded* cm er;

% Keep only agglomerates (including single segment agglomerados) over minSize voxel
axonsSize = cellfun(@(x)sum(segmentMeta.voxelCount(x)), axons);
axons(axonsSize < options.minSize) = [];

idx = 1:10000:length(axons);
axons = axons(idx);
result = connectEM.calculateDirectionalityOfAgglomerates(axons, graph, segmentMeta, borderMeta, globalSegmentPCA, options);

% Write some skeletons for quick check
for i=1:length(result.neighbours)
    clear comments;
    comments{1} = [];
    filename = ['/gaba/scratch/mberning/directionalityDebug/' num2str(i, '%.3i') '.nml'];
    theseNodes = cellfun(@(x)segmentMeta.point(:,x)', cat(1, axons(i), result.neighbours(i)), 'uni', 0);
    theseComments = arrayfun(@(x,y)[num2str(x, '%+3.2f') '  ' num2str(y, '%+3.2f')], result.latent{i}(:,1), result.scores{i}, 'uni', 0);  
    idx = result.latent{i}(:,1) > 0.7 & abs(result.scores{i}(:,1)) > 0.7;
    if sum(idx)
        comments{2} = theseComments(idx);
        nodes{1} = theseNodes{1};
        nodes{2} = theseNodes{2}(idx,:);
        connectEM.generateSkeletonFromNodes(filename, nodes, {'agglo' 'neighbours'}, comments, true);
    else
        connectEM.generateSkeletonFromNodes(filename, nodes(1), {'agglo'}, [], true);
    end
end

