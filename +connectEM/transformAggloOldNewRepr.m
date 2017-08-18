function aggloNew = transformAggloOldNewRepr(aggloOld,edgesSegId)

searchVec = cat(1,aggloOld{:});
[~,idx] = ismember(edgesSegId,searchVec);

numSegs = cellfun(@numel,aggloOld);

countVec = cat(1,NaN,repelem((1:numel(numSegs))',numSegs));  % create a indices vector to reference the idx to the agglo id in class1
aggloIdx = countVec(idx+1);  % make the agglomeration indices matrix corresponding to the corresponding edges

% check if there are no overlaps between agglos
assert(all(aggloIdx(:,1)==aggloIdx(:,2)))

% sort hot edges according to the agglos they will be put in
[saggloIdx,sIdx] = sort(aggloIdx(:,1));
edgesSegId = edgesSegId(sIdx,:);
% delete the edges not used (should not be the case)
edgesSegId(isnan(saggloIdx),:) = [];
saggloIdx(isnan(saggloIdx)) = [];

% transform the edges into cell and create node cell array including segID
% information
edges = mat2cell(edgesSegId,histc(saggloIdx,unique(saggloIdx)));
nodes = cellfun(@(x) [segmentMeta.point(:,x)' x],aggloOld,'uni',0);
% tranform the segIds in the edge vector to index to the node
[~, edges] = cellfun(@(x,y) ismember(x,y(:,4)),edges,nodes,'uni',0);

% transform into struct
aggloNew = cell2struct([edges';nodes'],{'edges','nodes'},1);

