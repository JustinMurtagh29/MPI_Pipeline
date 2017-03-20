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
graph = load([p.saveFolder 'graph.mat']);
% Load information about edges
borderMeta = load([p.saveFolder 'globalBorder.mat']);
% Load meta information of segments
segmentMeta = load([p.saveFolder 'segmentMeta.mat']);
%% Load synapse scores from SynEM
synScore = load([p.saveFolder 'globalSynScores.mat']);
synScore.isSynapse = connectEM.synScoresToSynEdges(graph, synScore);
toc;

display('Removing segments detected by heuristics:');
tic;
% Segments classified by heuristics
load([p.saveFolder 'heuristicResult.mat']);
assert(length(segIds) == max(segIds));
excludedIds1 = vesselScore > 0.5 | myelinScore > 0.5 | nucleiScore > 0.5;
% Keep only segments larger than 100 voxel
excludedIds2 = segmentMeta.voxelCount <= 100;
% Remove from graph
keepIdx = ~any(ismember(graph.edges, find(excludedIds1 | excludedIds2)),2);
graphCut.edges = graph.edges(keepIdx,:);
graphCut.prob = graph.prob(keepIdx);
toc;

rng default;
threshold=0.99:-0.01:0.91;
for t=1:5%:length(threshold)
    display('Performing agglomeration:');
    tic;
    %% Agglomerate segments using only GP probabilties
    [initialPartition{t}, remainingEdges] = connectEM.partitionWholeDataset(graphCut, threshold(t)); 
    sizePartition{t} = cellfun(@(x)sum(segmentMeta.voxelCount(x)), initialPartition{t});
    [sizePartition{t}, idx] = sort(sizePartition{t}, 'descend');
    initialPartition{t} = initialPartition{t}(idx);
    toc;
    % Probably make own function if needed again
    display('Writing skeletons for debugging the process:');
    tic;
    connectEM.generateSkeletonFromAgglo(remainingEdges, segmentMeta.point', ...
        initialPartition{t}(1:1000), ...
        strseq(['largestComponents' num2str(threshold(t)) '_'], 1:1000), ...
        outputFolder, segmentMeta.maxSegId);
    idx = randi(numel(initialPartition{t}),100,1);
    connectEM.generateSkeletonFromAgglo(remainingEdges, segmentMeta.point', ...
        initialPartition{t}(idx), ...
        strseq(['randomComponents' num2str(threshold(t)) '_'], 1:100), ...
        outputFolder, segmentMeta.maxSegId);
    toc;
end
voxelCount = segmentMeta.voxelCount;
maxSegId = segmentMeta.maxSegId;
save([outputFolder 'initialAgglo.mat'], 'initialPartition', 'sizePartition', 'voxelCount', 'maxSegId');


% FOR LATER
%{
%% Read query results
skeletonSaveFile = [outputFolder 'querySkeletonReadout.mat'];
if rereadSkeletons
    % Where to find skeletons that were returned from the queries
    scratchFolder = '/gaba/scratch/mberning/focusFlightTasks/extracted/';
    % These are the downloaded queries (see wK projects: L4_focus_flight-1 & 2, L4_focus_flight-reseed-1 & 2 (2 contains 2nd and 3rd round of reseeding now))
    % And new queries from after switching to new (agglomerate based) query analysis: L4_focus_flight-new-1
    skeletonFolders = {'L4_focus_flight-1', 'L4_focus_flight-2', 'L4_focus_flight-reseed-1', 'L4_focus_flight-reseed-2', ...
        'L4_focus_flight_new_1', 'L4_focus_flight_new_2' 'L4_focus_flight-new-3'};
    skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
    % Lookup segment ids of nodes+neighbours of nmls in all folders defined above
    [ff.segIds, ff.neighbours, ff.filenames, ff.nodes, ff.startNode] = lookupNmlMulti(p, skeletonFolders);
    % Read proofread subset (removed wrong seeds and proofread redundancies) of 100 5-fold subset
    [gt.segIds, gt.neighbours, gt.filenames, gt.nodes] = lookupNmlMulti(p, {'/gaba/scratch/mberning/100axons/training/'}, false);
    % Save for future reloading (above steps take hours otherwise (segmentation lookup)
    save(skeletonSaveFile, 'ff', 'gt');
else
    % Load skeleton readout including segmentation ids at nodes & all neighbours
    load(skeletonSaveFile, 'ff', 'gt');
end
toc;

display('Processing queries:');
tic;
%% Calculate overlap of all queries with segments 
[uniqueSegments, neighboursStartNode, nodesExcludedIdx, startNodeIdx] = cellfun(@queryAnalysis, ...
    ff.segIds, ff.neighbours, ff.nodes', ff.startNode', 'uni', 0); 
%% Determine all overlaps of agglomerations with given queries
[partition, queryOverlap] = queryAgglomerationOverlap(initialAgglomeration, uniqueSegments, neighboursStartNode);
%% Make decision(s), here evidence/occurence threshold is applied
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
%% Find CC of eqClasses to be joined
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

display('Generate new queries:');
tic;
% Decide which positions to querry and calculate some statistics
[q.pos, q.dir, q.varianceExplained, q.voxelSize, q.lengthAlongPC1] = cellfun(@(x)determineQueryLocation(graph, segments, x), initialAgglomeration.partition, 'uni', 0);
% Keep only axons which have a length along PC 1 direction >1 micron and 4000 voxel (40 voxel in z a 28 ~1 micron)
% Sort out all queries that are within 2 micron (query distance from border)
extend = round(2000 ./ [11.24 11.24 28]);
bbox(:,1) = p.bbox(:,1) + extend';
bbox(:,2) = p.bbox(:,2) - extend';
q.outsideBBox = cellfun(@(x)any(bsxfun(@le, x, bbox(:,1)'),2) | any(bsxfun(@ge, x, bbox(:,2)'),2), q.pos, 'uni', 0);
q.tooSmall = cellfun(@(x)x < 5e3, q.voxelSize);
q.tooShort = cellfun(@(x)x < 1.5e3, q.lengthAlongPC1);
q.tooUnstraight = cellfun(@(x)x(1) < .9, q.varianceExplained) & cellfun(@numel, initialAgglomeration.partition) > 2;
%nrBoutons = cellfun(@(x)sum(ismember(x, boutons.boutonIDs)), initialAgglomeration.partition);
%q.noSynapses = nrBoutons == 0;
q.exclude = arrayfun(@(x,y,z,t,u) x{:} | y | z | t, q.outsideBBox, q.tooShort, q.tooSmall, q.tooUnstraight, 'uni', 0);
% 'Write' problems to flight mode webKNOSSOS
extend = round(1000 ./ [11.24 11.24 28]);
fid = fopen([outputFolder 'webKnossos.txt'], 'w');
for i=1:length(q.pos)
    for j=1:size(q.pos{i},1)
        if ~q.exclude{i}(j)
            [phi, theta, psi] = calculateEulerAngles(-q.dir{i}(j,:));
            minPos = q.pos{i}(j,:) - extend;
            sizeBbox = 2*extend;
            linkString = ['2012-09-28_ex145_07x2,56d6a7c6140000d81030701e,focus_flight,1,' ...
                num2str(q.pos{i}(j,1)-1) ',' num2str(q.pos{i}(j,2)-1) ',' num2str(q.pos{i}(j,3)-1) ',' ...
                num2str(phi) ',' num2str(theta) ',' num2str(psi) ',1,Tracing crew,' ...
                num2str(minPos(1)) ',' num2str(minPos(2)) ',' num2str(minPos(3)) ',' ...
                num2str(sizeBbox(1)) ',' num2str(sizeBbox(2)) ',' num2str(sizeBbox(3)) ',' 'L4_focus_flight-new-4'];
            fprintf(fid, '%s\n', linkString);
            q.angles{i}(j,:) = [phi theta psi];
        end
    end
end
fclose(fid);
toc;

display('Writing skeletons for debugging the process:');
tic;
% TODO: Fix and make silent
% Write a debug & visualization for the queries overlaps (now with agglomerations instead of segments
debugQueryAttachment(graph, segments.pointMap, initialAgglomeration, ff, outputFolder);
% Write skeletons for comparing agglomerations merged by querries vs. ground truth
%debugAgglomerationToGtOverlap(graph, segments, gt, ff, partition, initialAgglomeration, boutons, eqClassCC, startAgglo, endAgglo, outputFolder);
debugAgglomerationToGtOverlap2(graph, segments, gt, initialAgglomeration, boutons, initialAgglomeration.partition, outputFolder);
debugAgglomerationToGtOverlap2(graph, segments, gt, initialAgglomeration, boutons, possibleAxons, outputFolder);
debugAgglomerationToGtOverlap3(graph, segments, gt, initialAgglomeration, boutons, partition, ff, startAgglo, endAgglo, outputFolder);
% Write skeletons for debugging new queries generated here in combination with agglomeration they are generated from
debugAgglomerationWithNewQueries(graph, segments, possibleAxons, q, outputFolder);
toc;

display('Saving all data in one humoungous file for bookkeeping purposes');
tic;
% Currently without flag which will loose us the graph variable (too large), but ok
% Add -v7.3 flag in case variables get to large?
save([outputFolder 'allData.mat']);
toc;

display('Calculating some statistics');
tic;
% How many (possible) axons with a certain number of boutons
for i=1:5
    display(['Number of axons with more than ' num2str(i) ' synapses: ' num2str(sum(nrBoutons > i))]);
end
nrEndsAtBorder = cellfun(@sum, q.outsideBBox);
lengthPC1 = cellfun(@(x)norm((x(1,:)-x(2,:)) .* [11.24 11.24 28]), q.pos);
for i=0:1
    display(['Number of axons with more than 10 micron path length and ' num2str(i) ' ends at dataset border: ' num2str(sum(nrEndsAtBorder > i & lengthPC1 > 10e3))]);
end
toc;
%}
