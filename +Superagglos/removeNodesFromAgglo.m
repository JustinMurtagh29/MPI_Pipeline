function [newSuperagglo] = removeNodesFromAgglo(superagglo,nodeIds )
% removes a set of segment IDs from all superagglos and performs connected
% component afterswards to ensure internode connectivity in each superagglo

nodes = cell2mat(arrayfun(@(x) x.nodes,superagglo,'uni',0));
fNames = setdiff(fieldnames(superagglo),{'nodes','edges'});
tmpstrct = cell(numel(fNames),1);
for f = 1:numel(fNames)
    tmpstrct{f} = cat(1,superagglo.(fNames{f}));
end
numNodes = arrayfun(@(x) size(x.nodes,1),superagglo);
cumNumNodes = [0;cumsum(numNodes)];
edges = cell2mat(arrayfun(@(x) superagglo(x).edges+cumNumNodes(x),(1:numel(superagglo))','uni',0));% putting all edges together (inceasing index for each superagglo

idsToKeep = setdiff((1:size(nodes,1))',nodeIds);
edges = edges(all(ismember(edges,idsToKeep),2),:);

% generate selfEdges of all superagglo segments
% selfEdges = repmat(nodes(:,4),1,2);
selfEdges = repmat(idsToKeep,1,2);

% [equivalenceClass, aggloLUT] = Graph.findConnectedComponents(cat(1,selfEdges,segIDedges),0,1);
[equivalenceClass, aggloLUT] = Graph.findConnectedComponents(cat(1,selfEdges,edges),0,1);
if isempty(equivalenceClass)
    newSuperagglo = struct('edges',[],'nodes',[]);
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
newedges = cell(numel(equivalenceClass),1);
% Treat agglomerates containing single segment separately
% newedges(~singleSegAgglos) = mat2cell(segIDedges,histc(saggloLUT,unique(saggloLUT)));
newedges(~singleSegAgglos) = mat2cell(edges,histc(saggloLUT,unique(saggloLUT)));
newedges(singleSegAgglos) = {zeros(0,2)}; % give empty edges correct dimension

%  create node cell array including segID information
newnodes = cellfun(@(x) nodes(x,:),equivalenceClass,'uni',0)';
for f = 1:numel(fNames)
    if size(tmpstrct{f},1) == size(nodes,1) % apply same sorting to all other field names if their size was the same as the number of nodes, else leave the same
        tmpstrct{f} = cellfun(@(x) tmpstrct{f}(x,:),equivalenceClass,'uni',0)';
    else
        tmpstrct{f} = repmat(tmpstrct(f),1,numel(equivalenceClass));
    end
end
% tranform the global node ids in the edge vector to local node ids
[~, newedges] = cellfun(@(x,y) ismember(x,y),newedges,equivalenceClass,'uni',0);

newSuperagglo = cell2struct([newedges;newnodes;tmpstrct{:}],[{'edges'},{'nodes'},fNames'],1);
end

