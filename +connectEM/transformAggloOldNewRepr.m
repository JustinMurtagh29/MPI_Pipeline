function aggloNew = transformAggloOldNewRepr(aggloOld,edgesSegId,segmentMeta)

searchVec = cat(1,aggloOld{:});
assert(all(ismember(searchVec,edgesSegId(:))))

numSegs = cellfun(@numel,aggloOld);
countVec = cat(1,NaN,repelem((1:numel(numSegs))',numSegs));  % create a indices vector to reference the idx to the agglo id in class1
while 1
    [~,idx] = ismember(edgesSegId,searchVec);
    aggloIdx = countVec(idx+1);  % make the agglomeration indices matrix corresponding to the corresponding edges
    
    % check if there are no overlaps between agglos
    if ~ all(aggloIdx(:,1)==aggloIdx(:,2))
        warning('Edge list also contained edges between agglos. These are thrown out now')
%         aggloIdx = aggloIdx(aggloIdx(:,1)==aggloIdx(:,2),:);
        edgesSegId(aggloIdx(:,1)~=aggloIdx(:,2),:) = [];
    else
        break
    end
end
% sort hot edges according to the agglos they will be put in
[saggloIdx,sIdx] = sort(aggloIdx(:,1));
edgesSegId = edgesSegId(sIdx,:);
% delete the edges not used (should not be the case)
edgesSegId(isnan(saggloIdx),:) = [];
saggloIdx(isnan(saggloIdx)) = [];

% transform the edges into cell and create node cell array including segID
% information
edges = cell(numel(numSegs),1);
edges(numSegs>1) = mat2cell(edgesSegId,histc(saggloIdx,unique(saggloIdx)));
edges(numSegs==1) = zeros(0,2); % give empty edges correct dimension
nodes = cellfun(@(x) [segmentMeta.point(:,x)' x],aggloOld,'uni',0);
% tranform the segIds in the edge vector to index to the node
[~, edges] = cellfun(@(x,y) ismember(x,y(:,4)),edges,nodes,'uni',0);

% transform into struct
aggloNew = cell2struct([edges';nodes'],{'edges','nodes'},1);

