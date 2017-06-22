
outputFolder = '/gaba/scratch/mberning/dendriteQueryResults/';

% Load data
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
segmentMeta = load([p.saveFolder 'segmentMeta.mat']);
segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);
temp = load('/gaba/scratch/mberning/aggloGridSearch/search03_00514.mat');
dendrites = temp.dendritesFinal;
clear temp
segmentsLeftover = setdiff(find(segmentMeta.dendriteProb > 0.5), cell2mat(dendrites));
% Not needed so far
%graph = load([p.saveFolder 'graphNew.mat'], 'edges', 'prob', 'borderIdx');

% Where to find skeletons that were returned from the queries
scratchFolder = '/gaba/scratch/mberning/queryResults2017/';
% These are the downloaded queries (see wK projects: L4_focus_flight-1 & 2, L4_focus_flight-reseed-1 & 2 (2 contains 2nd and 3rd round of reseeding now))
% And new queries from after switching to new (agglomerate based) query analysis: L4_focus_flight-new-1
skeletonFolders = {'MBKMB_L4_dendrite_queries_2017_a'};
skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
% Lookup segment ids of nodes+neighbours of nmls in all folders defined above
[ff.segIds, ff.neighbours, ff.filenames, ff.nodes, ff.startNode, ff.comments] = connectEM.lookupNmlMulti(p, skeletonFolders, false);
display([num2str(sum(~cellfun(@isempty,ff.comments))) '/' num2str(numel(ff.comments)) ' queries contain comment and will not be used']);
tabulate(cellfun(@(x)x{1}{1}, cellfun(@(x)regexp(x, 'content="(.*)"', 'tokens'), ...
    cat(1, ff.comments{~cellfun(@isempty, ff.comments)}), 'uni', 0), 'uni', 0))
% Where to find skeletons that were returned from the queries
scratchFolder = '/gaba/scratch/mberning/queryResults2017/';
ff = structfun(@(x)x(cellfun(@isempty,ff.comments)), ff, 'uni', 0);

% Calculate overlap of all queries with segments
[uniqueSegments, neighboursStartNode, nodesExcludedIdx, startNodeIdx] = cellfun(@connectEM.queryAnalysis, ...
        ff.segIds, ff.neighbours, ff.nodes', ff.startNode', 'uni', 0);
% Determine all overlaps of agglomerations with given queries
[partition, queryOverlap] = connectEM.queryAgglomerationOverlap(dendrites, segmentsLeftover, uniqueSegments, neighboursStartNode);
% Make decision(s), here evidence/occurence threshold is applied
% Always one (or none if evidence below 14, 1/2 node) start eqClass
startAgglo = arrayfun(@(x)x.eqClasses(x.occurences > 13), queryOverlap.start, 'uni', 0);
% Exclude all queries that do not have a clear starting point
idxNoClearStart = cellfun('isempty', startAgglo);
% Multiple ends (all above 53vx evidence, corresponds to 2 full nodes)
endAgglo = arrayfun(@(x)x.eqClasses(x.occurences > 53), queryOverlap.ends, 'uni', 0);
% Exclude startAgglo from endAgglo (as we do not want to count self-attachment)
endAgglo = cellfun(@(x,y)setdiff(x,y), endAgglo, startAgglo, 'uni', 0);
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

% Visualization of queries and connections made
idx = find(idxGood);
idx = idx(randperm(numel(idx), 50));
temp = structfun(@(x)x(idx), ff, 'uni', 0);
connectEM.debugQueryAttachment(segmentMeta.point', dendrites, temp, outputFolder, 'queryAttachted');
idx = find(~idxGood);
idx = idx(randperm(numel(idx), 50));
temp = structfun(@(x)x(idx), ff, 'uni', 0);
connectEM.debugQueryAttachment(segmentMeta.point', dendrites, temp, outputFolder, 'queryOpenEnd');
clear idx temp;

