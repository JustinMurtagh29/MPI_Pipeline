
%{
% load data
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat', 'voxelCount', 'centroid', 'box', 'maxSegId', 'cubeIdx', 'point');
temp = load('/gaba/scratch/kboerg/directCycle/cycle006.mat', 'axonsFinal', 'metrics');
agglos = temp.axonsFinal;
aggloSize = cellfun(@(x)sum(segmentMeta.voxelCount(x)), agglos);
[aggloSize, aggloIdx] = sort(aggloSize, 'descend');
agglos = agglos(aggloIdx);
agglosGtAxon2 = temp.axonsFinal(temp.metrics.axon1.foundAgglomerates_col{2});
agglosGtAllAxons = temp.axonsFinal(cat(1,temp.metrics.axon1.foundAgglomerates_col{1:10}));
[~, cm] = connectEM.getERcomponents();
excludedCubeIdx = unique(cellfun(@(x)mode(segmentMeta.cubeIdx(x)), cm));
excludedSegmentIdx = ismember(segmentMeta.cubeIdx, excludedCubeIdx);
% Single segments if not excluded due to catastrohic merger
segmentsLeftover = setdiff(find(segmentMeta.axonProb > 0.5 & ~excludedSegmentIdx), cell2mat(agglos));
clear excluded* cm temp;
% calculate queries and write to file on wk
outputFolder = '/gaba/scratch/mberning/flightQueries/';

% 2 different versions for MH (1 and 2 micron bounding box)
options.writeTasksToFile = true;
options.queryBoundingBoxSize = 1000;
options.datasetBorderExclusionSize = 2000;
q = connectEM.generateQueriesFromAgglos(p, segmentMeta, agglosGtAxon2, outputFolder, options);
connectEM.debugNewQueries(segmentMeta, agglosGtAxon2, q, outputFolder);
% 2nd version
options.writeTasksToFile = true;
options.queryBoundingBoxSize = 2000;
options.datasetBorderExclusionSize = 2000;
q = connectEM.generateQueriesFromAgglos(p, segmentMeta, agglosGtAxon2, outputFolder, options);
connectEM.debugNewQueries(segmentMeta, agglosGtAxon2, q, outputFolder);

% Do it on whole dataset for query number vs. agglomerateSize plots
options.writeTasksToFile = false;
options.queryBoundingBoxSize = 1000;
options.datasetBorderExclusionSize = 2000;
q = connectEM.generateQueriesFromAgglos(p, segmentMeta, agglos, outputFolder, options);
save([outputFolder 'allQueries.mat'], 'q');

% All GT axon for query analysis tests
options.writeTasksToFile = true;
options.queryBoundingBoxSize = 2000;
options.datasetBorderExclusionSize = 2000;
q = connectEM.generateQueriesFromAgglos(p, segmentMeta, agglosGtAllAxons, outputFolder, options);
connectEM.debugNewQueries(segmentMeta, agglosGtAllAxons, q, outputFolder);
%}

%% Query results first tests
% Where to find skeletons that were returned from the queries
scratchFolder = '/gaba/scratch/mberning/';
% These are the downloaded queries (see wK projects: L4_focus_flight-1 & 2, L4_focus_flight-reseed-1 & 2 (2 contains 2nd and 3rd round of reseeding now))
% And new queries from after switching to new (agglomerate based) query analysis: L4_focus_flight-new-1
skeletonFolders = {'axonQueries'};
skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
% Lookup segment ids of nodes+neighbours of nmls in all folders defined above
[ff.segIds, ff.neighbours, ff.filenames, ff.nodes, ff.startNode, ff.comments] = connectEM.lookupNmlMulti(p, skeletonFolders);
display([num2str(sum(~cellfun(@isempty,ff.comments))) '/' num2str(numel(ff.comments)) ' queries contain comment and will not be used']);
ff = structfun(@(x)x(cellfun(@isempty,ff.comments)), ff, 'uni', 0);

% Intermediate visualization
connectEM.debugQueryAttachment(segmentMeta.point', agglos, ff, outputFolder);

display('Processing queries:');
tic;
% Calculate overlap of all queries with segments
[uniqueSegments, neighboursStartNode, nodesExcludedIdx, startNodeIdx] = cellfun(@connectEM.queryAnalysis, ...
    ff.segIds, ff.neighbours, ff.nodes', ff.startNode', 'uni', 0); 
% Determine all overlaps of agglomerations with given queries
[partition, queryOverlap] = connectEM.queryAgglomerationOverlap(agglos, segmentsLeftover, uniqueSegments, neighboursStartNode);
% Make decision(s), here evidence/occurence threshold is applied
% Always one (or none if evidence below 14, 1/2 node) start eqClass
startAgglo = arrayfun(@(x)x.eqClasses(x.occurences > 13), queryOverlap.start, 'uni', 0);
% Exclude all queries that do not have a clear starting point (6.5% on first 200k 07x2)
idxNoClearStart = cellfun('isempty', startAgglo);
% Multiple ends (all above 53vx evidence, corresponds to 2 full nodes)
endAgglo = arrayfun(@(x)x.eqClasses(x.occurences > 53), queryOverlap.ends, 'uni', 0);
% Exclude all queries that do not have (at least one) clear end
idxNoClearEnd = cellfun('isempty', endAgglo); 
% 18.5% of queries excluded overall due to missing start or end (or both)
idxGood = ~(idxNoClearStart | idxNoClearEnd);
% Display some statistics
display([num2str(sum(idxNoClearStart)./numel(idxNoClearStart)*100, '%.2f') '% of remaining queries have no clear start']);
display([num2str(sum(idxNoClearEnd)./numel(idxNoClearEnd)*100, '%.2f') '% of remaining queries have no clear end']);
display([num2str(sum(idxGood)./numel(idxGood)*100, '%.2f') '% of remaining queries have clear start and ending']);
display([num2str(numel(cat(2, endAgglo{idxGood}))) ' attachments made by ' num2str(sum(idxGood)) ' queries']);
% Find CC of eqClasses to be joined
edges = cellfun(@(x,y)combnk([x y], 2), startAgglo(idxGood), endAgglo(idxGood), 'uni', 0);
edges = cat(1,edges{:});
edges(edges(:,1) == edges(:,2),:) = [];
eqClassCC = Graph.findConnectedComponents(edges, true, true);
toc;

display('Merging agglomerates based on queries:');
tic;
% Build data structure for faster (?) merging and merge
partitionNew.segIds = partition.segIds;
partitionNew.eqClass = partition.eqClass;
for i=1:length(eqClassCC)
    partitionNew.eqClass(ismember(partitionNew.eqClass, eqClassCC{i}(2:end))) = eqClassCC{i}(1);
end
[partitionNew.eqClass, idx] = sort(partitionNew.eqClass);
partitionNew.segIds = partitionNew.segIds(idx);
eqClassEnds = cat(1, find(diff(partitionNew.eqClass)), length(partitionNew.eqClass));
eqClassSize = cat(1, eqClassEnds(1), eqClassEnds(2:end) - eqClassEnds(1:end-1));
possibleAxons = mat2cell(partition.segIds, eqClassSize);
toc;


