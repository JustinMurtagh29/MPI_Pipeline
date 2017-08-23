function aggloNew = transformAggloOldNewRepr(aggloOld,edgesSegId,segmentMeta,doChecks)

if ~exist('doChecks','var')
    doChecks = 0;
end

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
edges(numSegs>1) = mat2cell(edgesSegId,histc(sedgesLookup,unique(sedgesLookup)));
edges(numSegs==1) = {zeros(0,2)}; % give empty edges correct dimension
nodes = cellfun(@(x) [segmentMeta.point(:,x)' x],aggloOld,'uni',0);
% tranform the segIds in the edge vector to index to the node
[~, edges] = cellfun(@(x,y) ismember(x,y(:,4)),edges,nodes,'uni',0);

% transform into struct
aggloNew = cell2struct([edges';nodes'],{'edges','nodes'},1);

if doChecks
    % Perform some sanity checks
    idxSingleSegment = numSegs == 1;
    
    % Check that each agglo is own CC
    assert(all(arrayfun(@(x)numel(Graph.findConnectedComponents(x.edges)), ...
        aggloNew(~idxSingleSegment)) == 1));
    
    % Check that all nodes are connected
    assert(all(arrayfun(@(x)isempty(setdiff(1:size(x.nodes,1),x.edges(:))), ...
        aggloNew(~idxSingleSegment))));
    assert(all(arrayfun(@(x)isempty(setdiff(x.edges(:),1:size(x.nodes,1))), ...
        aggloNew(~idxSingleSegment))));
    
    % Check that equivalence classes are exclusive
    segIds = cell2mat(arrayfun(@(x)x.nodes(:,4), aggloNew, 'uni', 0));
    assert(all(histc(segIds, unique(segIds)) == 1));
    
    % Check that segments still belong to the same equivalence class
    assert(all(arrayfun(@(x)isempty(setdiff(aggloOld{x},aggloNew(x).nodes(:,4))), 1:numel(aggloOld))));
    assert(all(arrayfun(@(x)isempty(setdiff(aggloNew(x).nodes(:,4),aggloOld{x})), 1:numel(aggloOld))));
end
end
