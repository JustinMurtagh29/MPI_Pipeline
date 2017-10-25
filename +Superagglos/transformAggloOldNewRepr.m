function aggloNew = transformAggloOldNewRepr(aggloOld,edgesSegId,segmentMeta,doChecks)

if ~exist('doChecks','var')
    doChecks = 0;
end

% make sort the edges in second dim and make them unique
edgesSegId = unique(sort(edgesSegId,2),'rows');

% Lookup table to which agglo each segment id belongs
numSegs = cellfun(@numel, aggloOld);
aggloLookup = NaN(segmentMeta.maxSegId,1);
aggloLookup(cell2mat(aggloOld)) = ...
    repelem(1:numel(aggloOld), numSegs);

% Lookup table which edge belongs to which agglos
edgesLookup = aggloLookup(edgesSegId);
intraIdx = edgesLookup(:,1) == edgesLookup(:,2);

% Filter out non-intra agglo edges
edgesSegId = edgesSegId(intraIdx,:);
edgesLookup = edgesLookup(intraIdx,1);
assert(~any(isnan(edgesLookup(:))));

% sort hot edges according to the agglos they will be put in
[sedgesLookup,sIdx] = sort(edgesLookup);
edgesSegId = edgesSegId(sIdx,:);

% transform the edges into cell and create node cell array including segID
% information
edges = cell(numel(aggloOld),1);
% Treat agglomerates containing single segment separately
try
    edges(numSegs>1) = mat2cell(edgesSegId,histc(sedgesLookup,unique(sedgesLookup)));
catch
    error('Error in edge assignment. Please check if an agglo has members that are not interconnected via the given edges')
end
edges(numSegs==1) = {zeros(0,2)}; % give empty edges correct dimension
nodes = cellfun(@(x) [segmentMeta.point(:,x)' x],aggloOld,'uni',0);
% tranform the segIds in the edge vector to index to the node
[~, edges] = cellfun(@(x,y) ismember(x,y(:,4)),edges,nodes,'uni',0);

% transform into struct
aggloNew = cell2struct([edges';nodes'],{'edges','nodes'},1);

if doChecks
    % Perform some sanity checks
    idxSingleSegment = numSegs == 1;
      
    % Check that all nodes are connected
    assert(all(arrayfun(@(x)isempty(setxor(1:size(x.nodes,1),x.edges(:))), ...
        aggloNew(~idxSingleSegment))));
    
    % Check that equivalence classes are exclusive
    segIds = cell2mat(arrayfun(@(x)x.nodes(:,4), aggloNew, 'uni', 0));
    assert(all(histc(segIds, unique(segIds)) == 1));
    
    % Check that segments still belong to the same equivalence class
    assert(all(arrayfun(@(x)isempty(setxor(aggloOld{x},aggloNew(x).nodes(:,4))), 1:numel(aggloOld))));
    
    % Check that there are no edge duplets in any agglo
    assert(all(arrayfun(@(x) size(reshape(x.nodes(x.edges,4),[],2),1)==size(unique(sort(reshape(x.nodes(x.edges,4),[],2),2),'rows'),1),aggloNew(~idxSingleSegment))))
end
end
