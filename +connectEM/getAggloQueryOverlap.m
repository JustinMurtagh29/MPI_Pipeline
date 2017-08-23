function getAggloQueryOverlap(superagglos)

load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');

scratchFolder = '/tmpscratch/mberning/axonQueryResults/';
skeletonFolders = {'MBKMB_L4_axons_queries_2017_a' 'MBKMB_L4_axons_queries_2017_b'};
skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
skeletonFolders = [skeletonFolders {'/tmpscratch/scchr/AxonEndings/axonQueryResults/CS_MB_L4_AxonLeftQueries_nmls/'}];

for i=1:length(superagglos)
    axons{i,1} = superagglos(i).nodes(:,4);
end
axons = cellfun(@(x)x(~isnan(x)),axons,'uni',0);

% Lookup segment ids of nodes+neighbours of nmls in all folders defined above
[ff.segIds, ff.neighbours, ff.filenames, ff.nodes, ff.startNode, ff.comments] = connectEM.lookupNmlMulti(p, skeletonFolders, false);

display([num2str(sum(~cellfun(@isempty,ff.comments))) '/' num2str(numel(ff.comments)) ' queries contain comment and will not be used']);
tabulate(cellfun(@(x)x{1}{1}, cellfun(@(x)regexp(x, 'content="(.*)"', 'tokens'), ...
    cat(1, ff.comments{~cellfun(@isempty, ff.comments)}), 'uni', 0), 'uni', 0))

% queries commented by the HIWI in terms of unsolvability
idx_comment = ~cellfun(@isempty,ff.comments);

% ~600 queries do not have a start node, not sure why (maybe the ones with more than one tree), maybe check later
idx_startEmpty = cellfun(@isempty, ff.startNode);

save([p.saveFolder 'AxonFlightPaths.mat'], 'ff', 'idx_comment', 'idx_startEmpty')
% load([p.saveFolder 'AxonFlightPaths.mat'])

ff = structfun(@(x)x(cellfun(@isempty, ff.comments)), ff, 'uni', 0);
ff = structfun(@(x)x(~cellfun(@isempty, ff.startNode)), ff, 'uni', 0);
%%    
segmentsLeftover = [];

% Calculate overlap of all queries with segments
[uniqueSegments, neighboursStartNode, nodesExcludedIdx, startNodeIdx] = cellfun(@connectEM.queryAnalysis, ...
        ff.segIds, ff.neighbours, ff.nodes', ff.startNode', 'uni', 0);
% Determine all overlaps of agglomerations with given queries
[partition, queryOverlap] = connectEM.queryAgglomerationOverlap(axons, segmentsLeftover, uniqueSegments, neighboursStartNode);
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
% Find CC of eqClasses to be joined including single eqClasses with or
% without dangling query
edgesCC = cellfun(@(x,y)combnk([x y], 2), startAgglo(idxGood), endAgglo(idxGood), 'uni', 0);
edgesCC = cat(1,edgesCC{:});
edgesCC = sort(edgesCC, 2);
eqClassCC = Graph.findConnectedComponents(edgesCC, true, true);
sizeEqClassCC = sort(cellfun(@numel, eqClassCC), 'descend');
eqClassCCfull = [eqClassCC; num2cell(setdiff(1 : length(axons), cell2mat(eqClassCC)))'];
display(sizeEqClassCC(1:10));

results.startAgglo = startAgglo;
results.endAgglo = endAgglo;
results.ff = ff;
results.idxGood = idxGood;

save([p.savefolder 'AxonQueryOverlaps.mat'], 'results', 'queryOverlap', 'idxNoClearStart', 'idxNoClearEnd')
save([p.savefolder 'AxonPostQueryAnalysisState.mat'])
