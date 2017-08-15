outputFolder = '/tmpscratch/mberning/axonQueryResults/';

load('/gaba/scratch/mberning/edgesGTall.mat');
% Load segment meta data & graph
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
[graph, segmentMeta] = connectEM.loadAllSegmentationData(p);
% Load state of axons > 5 micron
temp = load('/gaba/scratch/mberning/axonQueryGeneration/beforeQueryGeneration.mat', 'axonsNew');
axons = temp.axonsNew;
clear temp;
functionH = @connectEM.axonQueryAnalysisSub;
inputCell = cellfun(@(x){x}, num2cell(1 : 500), 'uni', 0);
cluster = Cluster.getCluster( ...
    '-pe openmp 1', ...
    '-p 0', ...
    '-l h_vmem=24G', ...
    '-l s_rt=23:50:00', ...
    '-l h_rt=24:00:00');
job = Cluster.startJob( functionH, inputCell, ...
    'name', 'query', ...
    'cluster', cluster);
for startidx = 1  : 500
    startidx
    temp = load(['/tmpscratch/kboerg/20170810axonQueryAnalysis/output' num2str(startidx)],'usededges');
    usededges(startidx:500:size(axons,1)) = temp.usededges(startidx:500:size(axons,1));
end
for idx = 1 : length(axons)
    idx
    assert(isequal(Graph.findConnectedComponents(edgesGTall(usededges{idx},:)),{sort(axons{idx})}));
end
usededges2=cell2mat(usededges(~cellfun(@isempty,usededges))');
usededges3 =ismember(graph.edges,edgesGTall(usededges2,:),'rows');
% Only use 5 micron pieces in this initial analysis
%segmentsLeftover = setdiff(find(segmentMeta.axonProb > 0.5), cell2mat(axons));
segmentsLeftover = [];

options.border = [2000; -2000];
borderNm = repmat(options.border, 1, 3);
borderVoxel = round(bsxfun(@times, 1./p.raw.voxelSize, borderNm));
bboxSmall = p.bbox + borderVoxel';
borderMeta = load([p.saveFolder 'globalBorder.mat'], 'borderCoM');
borderMeta.borderCoM(end+1,:)=nan(1,3);
borderidxs = graph.borderIdx(usededges3);
borderidxs(isnan(borderidxs))=length(borderMeta.borderCoM);

borderPositions = double(borderMeta.borderCoM(borderidxs,:));

outsideBbox = ~(all(bsxfun(@gt, borderPositions, bboxSmall(:, 1)'), 2) & ...
    all(bsxfun(@lt, borderPositions, bboxSmall(:, 2)'), 2));
usededges4=usededges3;    
usededges4(graph.prob(usededges3)<0.98&outsideBbox)  = false;
axonsNew = Graph.findConnectedComponents(graph.edges(usededges4,:));
axonsNew = [axonsNew; num2cell(setdiff(cell2mat(axons), cell2mat(axonsNew)))];
calculateLength = @(x)max(pdist(bsxfun(@times, double(borderMeta.borderCoM(x, :)), p.raw.voxelSize)));
filternan = @(x)x(~isnan(x));
for idx = 1 : length(axonsNew)
    idx
    axonLength(idx) = max([-1, calculateLength(filternan(cell2mat(graph.neighBorderIdx(axonsNew{idx}))))]);
end
axons=axonsNew(axonLength>5000);
% Where to find skeletons that were returned from the queries
scratchFolder = outputFolder;
skeletonFolders = {'MBKMB_L4_axons_queries_2017_a' 'MBKMB_L4_axons_queries_2017_b'};
skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);
% Lookup segment ids of nodes+neighbours of nmls in all folders defined above
[ff.segIds, ff.neighbours, ff.filenames, ff.nodes, ff.startNode, ff.comments] = connectEM.lookupNmlMulti(p, skeletonFolders, false);
display([num2str(sum(~cellfun(@isempty,ff.comments))) '/' num2str(numel(ff.comments)) ' queries contain comment and will not be used']);
tabulate(cellfun(@(x)x{1}{1}, cellfun(@(x)regexp(x, 'content="(.*)"', 'tokens'), ...
    cat(1, ff.comments{~cellfun(@isempty, ff.comments)}), 'uni', 0), 'uni', 0))
ff = structfun(@(x)x(cellfun(@isempty,ff.comments)), ff, 'uni', 0);
% ~600 queries do not have a start node, not sure why (maybe the ones with more than one tree), maybe check later
ff = structfun(@(x)x(~cellfun(@isempty, ff.startNode)), ff, 'uni', 0);

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
% Find CC of eqClasses to be joined
edges = cellfun(@(x,y)combnk([x y], 2), startAgglo(idxGood), endAgglo(idxGood), 'uni', 0);
edges = cat(1,edges{:});
edges = sort(edges, 2);
eqClassCC = Graph.findConnectedComponents(edges, true, true);
sizeEqClassCC = sort(cellfun(@numel, eqClassCC), 'descend');
display(sizeEqClassCC(1:10));

% Generate new merged classes
eqClassCCfull = [eqClassCC; num2cell(setdiff(1 : length(axons), cell2mat(eqClassCC)))'];
axonsNew = cellfun(@(x){cell2mat(axons(x))}, eqClassCCfull);
% Load sorted out (and not used here) small agglomerates and add for Alessandro
load('/gaba/scratch/mberning/axonQueryGeneration/axonsSmall.mat', 'axonsSmall');
axonsPostQuery = cat(1, axonsNew, axonsSmall);
save([outputFolder 'postQueryAgglomerates.mat'], 'axonsPostQuery');

%Save workspace
save([outputFolder 'postQueryState.mat']);

%{
% Count unique rows
[unique_rows,~,ind] = unique(edges,'rows');
unique_ind = unique(ind);
counts = histc(ind,unique_ind);
tabulate(counts)
% Keep only edges that were voted for at least twice
edgesAbove1 = edges(unique_rows(unique_ind(counts > 1)),:);
eqClassCCabove1 = Graph.findConnectedComponents(edgesAbove1, true, true);
sizeEqClassCCAbove1 = sort(cellfun(@numel, eqClassCCabove1), 'descend');
display(sizeEqClassCCAbove1(1:10));

% Visualization of queries and connections made
idx = randperm(numel(ff.segIds), 100);
temp = structfun(@(x)x(idx), ff, 'uni', 0);
connectEM.debugQueryAttachment(segmentMeta.point', axons, temp, outputFolder, 'query');
clear idx temp;

% Generate skeletons for largest component for debuging
sizeEqClassCC = cellfun(@numel, eqClassCC);
[maxCC, maxCCidx] = max(sizeEqClassCC);
if ~exist([outputFolder 'percolator'], 'dir')
    mkdir([outputFolder 'percolator']);
end
skeletonFromAgglo(graph.edges, segmentMeta, axons(eqClassCC{maxCCidx}), 'mergedAgglos', [outputFolder 'percolator' filesep]);
%}

%{
% Look at redundant control tasks
liste = ff.filenames;
taskstridx=122;
taskStrings = unique(cellfun(@(x){x(1:taskstridx)}, liste))';
taskStrings{1,2} = [];
getUname = @(x)cellfun(@(y){y(taskstridx+1:end-12)}, x);
for idx = 1 : length(liste)
    if mod(idx, 100) == 0
        idx
    end
    rowidx = strcmp(taskStrings(:, 1), liste{idx}(1:taskstridx));
    taskStrings{rowidx, 2} = [taskStrings{rowidx, 2}, idx];
end
unames = unique(getUname(liste))';
unames{1, 2} = [];
testTasks = find(cellfun('length', taskStrings(:, 2)) == 3);
isequal2 = @(x,y)isequal(x,y)||(isempty(x)&&isempty(y));
allEqual = @(x)all(cellfun(@(y)isequal2(x{1}, y), x(2:end)));
for idx = cell2mat(taskStrings(testTasks, 2))'
    if allEqual(startAgglo(idx)) && allEqual(endAgglo(idx))
        uidxs = ismember(unames(:, 1), getUname(ff.filenames(idx)));
        unames(uidxs, 2) = cellfun(@(x){[x, 1]}, unames(uidxs, 2));
    else
        for badidx = idx'
            goodidx = setdiff(idx, badidx);
            if allEqual(startAgglo(goodidx)) && allEqual(endAgglo(goodidx))
                uidxs = ismember(unames(:, 1), getUname(ff.filenames(badidx)));
                unames(uidxs, 2) = cellfun(@(x){[x, -1]}, unames(uidxs, 2));
            end
        end
    end
end
%}
