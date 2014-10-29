% 
saveLocation = '/zdata/manuel/sync/activeTraining/';
load([saveLocation 'classData.mat']);

% make decision on error/reject of edges
errorCost = [1 5 10 50 100 500 1000 5000 10000 50000 100000]; % cost of error in relation to rejection
for i=1:length(errorCost)
    [segNew{i}, edgesNew{i}, pNew{i}] = joinSegments(seg, edges, p, errorCost(i));
    segNew{i} = imdilate(segNew{i}, ones(3,3,3));
end
save([saveLocation 'errorCostSeries.mat']);

%% Resort
errorCostIdx = 9;
[pResort, idx] = sort(pNew{errorCostIdx}, 'descend');
edgesResort = edgesNew(idx,:);

% transfer everything to knowledgeDB
segMerged(257:end-256,257:end-256,129:end-128) = segNew{errorCostIdx};
[settings, missions] = writeKnowledgeDB(paramBG, segMerged, edgesResort, pResort);
