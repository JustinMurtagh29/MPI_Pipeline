
outputFolder = '/gaba/scratch/kboerg/axonQueryResultsForHeiko/';

% Load data
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
segmentMeta = load([p.saveFolder 'segmentMeta.mat']);
segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);
temp = load('/gaba/scratch/mberning/aggloGridSearch6/6_01_00046/metricsFinal.mat', 'axonsNew');
axons = temp.axonsNew;
clear temp

% Not needed so far
%graph = load([p.saveFolder 'graphNew.mat'], 'edges', 'prob', 'borderIdx');

% Where to find skeletons that were returned from the queries
scratchFolder = '/gaba/scratch/kboerg/';
% These are the downloaded queries (see wK projects: L4_focus_flight-1 & 2, L4_focus_flight-reseed-1 & 2 (2 contains 2nd and 3rd round of reseeding now))
% And new queries from after switching to new (agglomerate based) query analysis: L4_focus_flight-new-1
skeletonFolders = {'annotationZips6791002912583435199HW_L4_GroundTruthForFlightTraining_nmls'};
skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0); 
% Lookup segment ids of nodes+neighbours of nmls in all folders defined above
[ff.segIds, ff.neighbours, ff.filenames, ff.nodes, ff.startNode, ff.comments] = connectEM.lookupNmlMulti(p, skeletonFolders, false);
ff.segIds(cellfun('isempty', ff.startNode)) = [];
ff.neighbours(cellfun('isempty', ff.startNode)) = [];
ff.filenames(cellfun('isempty', ff.startNode)) = [];
ff.nodes(cellfun('isempty', ff.startNode)) = [];
ff.comments(cellfun('isempty', ff.startNode)) = [];
ff.startNode(cellfun('isempty', ff.startNode)) = [];
display([num2str(sum(~cellfun(@isempty,ff.comments))) '/' num2str(numel(ff.comments)) ' queries contain comment and will not be used']);
tabulate(cellfun(@(x)x{1}{1}, cellfun(@(x)regexp(x, 'content="(.*)"', 'tokens'), ...
    cat(1, ff.comments{~cellfun(@isempty, ff.comments)}), 'uni', 0), 'uni', 0))
% Where to find skeletons that were returned from the queries
%scratchFolder = '/gaba/scratch/mberning/queryResults2017/';
ff = structfun(@(x)x(cellfun(@isempty,ff.comments)), ff, 'uni', 0); 

% Calculate overlap of all queries with segments
[uniqueSegments, neighboursStartNode, nodesExcludedIdx, startNodeIdx] = cellfun(@connectEM.queryAnalysis, ...
        ff.segIds, ff.neighbours, ff.nodes', ff.startNode', 'uni', 0);
% Determine all overlaps of agglomerations with given queries
[partition, queryOverlap] = connectEM.queryAgglomerationOverlap(axons, [], uniqueSegments, neighboursStartNode);
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
% idx = find(idxGood);
% idx = idx(randperm(numel(idx), 50));
% temp = structfun(@(x)x(idx), ff, 'uni', 0); 
% connectEM.debugQueryAttachment(segmentMeta.point', dendrites, temp, outputFolder, 'queryAttachted');
% idx = find(~idxGood);
% idx = idx(randperm(numel(idx), 50));
% temp = structfun(@(x)x(idx), ff, 'uni', 0); 
% connectEM.debugQueryAttachment(segmentMeta.point', dendrites, temp, outputFolder, 'queryOpenEnd');
% clear idx temp;
% 
% % Make complete list
% 
% eqClassCCfull = [eqClassCC; num2cell(setdiff(1 : length(axons), cell2mat(eqClassCC)))'];
% axonsPostQuery = cellfun(@(x){cell2mat(axons(x))}, eqClassCCfull);
% y_pre = connectEM.evaluateAggloMetaMeta(graph, axons, [],  'postQueryPreAxons', segmentMeta);
% y_post = connectEM.evaluateAggloMetaMeta(graph, axonsPostQuery,[],  'postQueryPostAxons', segmentMeta);
% 
% cell2mat(ff.startNode(find(cellfun(@(x)max([0,ismember(x,y_pre.axon1.foundAgglomerates_col{2})]), startAgglo)))')
% ff.filenames(find(cellfun(@(x)max([0,ismember(x,y_pre.axon1.foundAgglomerates_col{1})]), endAgglo)))'
% last=@(x)x{end};

% for idx = 1 : length(y_pre.axon1.foundAgglomerates_col)
%     folder = ['/gaba/scratch/kboerg/eval_agglo/postQueryPostAxons/axons1/query_' num2str(idx) '_'];
%     usedAgglos = cell2mat(eqClassCC(cellfun(@(x)any(ismember(x, y_pre.axon1.foundAgglomerates_col{idx})),eqClassCC)));
%     filenames = ff.filenames(cellfun(@(x,y)any(ismember([x,y],usedAgglos)), startAgglo, endAgglo));
%     
%     for idx2 = 1 : length(filenames)
%         copyfile(filenames{idx2}, [folder, last(strsplit(filenames{idx2}, filesep))]);
%     end
% end
liste = ff.filenames;
startidx = 148
taskStrings = unique(cellfun(@(x){x(1:startidx)}, liste))';
taskStrings{1,2} = [];
getUname = @(x)cellfun(@(y){y(startidx+1:end-10)}, x);
unames = unique(getUname(liste))';
unames{1, 2} = [];
for idx = 1 : length(liste)
    if mod(idx, 100) == 0
        idx
    end
    rowidx = strcmp(taskStrings(:, 1), liste{idx}(1:startidx));
    taskStrings{rowidx, 2} = [taskStrings{rowidx, 2}, idx];
end 
testTasks = find(cellfun('length', taskStrings(:, 2)) >= 3);
allEqual = @(x)all(cellfun(@(y)isequal(x{1}, y), x(2:end)));
memorizeTask=[];
for idx2 = 1:length(testTasks)
    idx = taskStrings{idx2, 2};
    if allEqual(startAgglo(idx)) && allEqual(endAgglo(idx))
        uidxs = ismember(unames(:, 1), getUname(ff.filenames(idx)));
        unames(uidxs, 2) = cellfun(@(x){[x, 1]}, unames(uidxs, 2));
        if length(idx) > 8
            memorizeTask(end+1)=idx2;
            
        end
    else
        for badidx = idx'
            goodidx = setdiff(idx, badidx);
            if allEqual(startAgglo(goodidx)) && allEqual(endAgglo(goodidx))
               uidxs = ismember(unames(:, 1), getUname(ff.filenames(badidx)));
               unames(uidxs, 2) = cellfun(@(x){[x, -1]}, unames(uidxs, 2));
               if length(idx) > 9
                   memorizeTask(end+1)=idx2;
               end
               break;
            end
        end
    end
end
unames(:, 2)= cellfun(@mat2str, unames(:, 2), 'uni', 0);
batch = load(['/gaba/scratch/mberning/axonQueryGenerationReverse/queriesMat/batch0050.mat']);
assert(length(unique(cell2mat(batch.q.pos),'rows' ))==507)
idx_batch=[];
for idx =1  : length(memorizeTask)
    idx_batch(end+1) = find(ismember(cell2mat(batch.q.pos)+1,ff.startNode{taskStrings{memorizeTask(idx),2}(1)},'rows'));
end
batch.q.pos=batch.q.pos(idx_batch);
batch.q.dir=batch.q.dir(idx_batch);
batch.q.angles=batch.q.angles(idx_batch);
save('/gaba/scratch/kboerg/HiwiBatch.mat','batch')
