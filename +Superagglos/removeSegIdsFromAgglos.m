function [newSuperagglos] = removeSegIdsFromAgglos(superagglos,ids )
% removes a set of segment IDs from all superagglos and performs connected
% component afterswards to ensure internode connectivity in each superagglo

nodes = cell2mat(arrayfun(@(x) x.nodes,superagglos,'uni',0));

numNodes = arrayfun(@(x) size(x.nodes,1),superagglos);
cumNumNodes = [0;cumsum(numNodes)];
edges = cell2mat(arrayfun(@(x) superagglos(x).edges+cumNumNodes(x),(1:numel(superagglos))','uni',0));% putting all edges together (inceasing index for each superagglo

% segIDedges = cell2mat(arrayfun(@(x) reshape(x.nodes(x.edges,4),[],2),superagglos,'uni',0));
idsToKeep = find(~ismember(nodes(:,4),ids));
edges = edges(all(ismember(edges,idsToKeep),2),:);

% nodeLUT(nodes(:,4)) = 1:size(nodes,1);
% segIDedges = segIDedges(all(~ismember(segIDedges,ids),2),:);

% generate selfEdges of all superagglo segments
% selfEdges = repmat(nodes(:,4),1,2);
selfEdges = repmat(idsToKeep,1,2);

% [equivalenceClass, aggloLUT] = Graph.findConnectedComponents(cat(1,selfEdges,segIDedges),0,1);
[equivalenceClass, aggloLUT] = Graph.findConnectedComponents(cat(1,selfEdges,edges),0,1);
if isempty(equivalenceClass)
    newSuperagglos = struct('edges',[],'nodes',[]);
    return
end
% create boolean which equivalence classes contain single edges
singleSegAgglos = true(numel(equivalenceClass),1);
% singleSegAgglos(aggloLUT(segIDedges(:))) = false;
singleSegAgglos(aggloLUT(edges(:))) = false;

% create lookup table which edge belongs to which
% equivalence class
% edgesLUT = aggloLUT(segIDedges(:,1));
edgesLUT = aggloLUT(edges(:,1));

% sort this LUT according to the agglos they will be put in and sort
% edges accordingly
[saggloLUT,sIdx] = sort(edgesLUT);
% segIDedges = segIDedges(sIdx,:);
edges = edges(sIdx,:);

% transform the edges into cell
newedges = cell(numel(superagglos),1);
% Treat agglomerates containing single segment separately
% newedges(~singleSegAgglos) = mat2cell(segIDedges,histc(saggloLUT,unique(saggloLUT)));
newedges(~singleSegAgglos) = mat2cell(edges,histc(saggloLUT,unique(saggloLUT)));
newedges(singleSegAgglos) = {zeros(0,2)}; % give empty edges correct dimension

%  create node cell array including segID information
newnodes = cellfun(@(x) nodes(x,:),equivalenceClass,'uni',0)';

% tranform the global node ids in the edge vector to local node ids
[~, newedges] = cellfun(@(x,y) ismember(x,y(:,4)),newedges,newnodes,'uni',0);

newSuperagglos = cell2struct([newedges;newnodes],{'edges','nodes'},1);
end

