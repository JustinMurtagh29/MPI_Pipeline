function [newSuperagglos] = removeSegIdsFromAgglos(superagglos,ids, noCC )
% removes a set of segment IDs from all superagglos and, if noCC flag is not true performs connected
% component afterwards to ensure internode connectivity in each superagglo
if ~exist('noCC','var') || isempty(noCC)
    noCC = 0;
end
nodes = cell2mat(arrayfun(@(x) x.nodes,superagglos,'uni',0));
toKeep = (~ismember(nodes(:,4),ids));
idsToKeep = find(toKeep);
idsToDelete = find(~toKeep);
numNodes = arrayfun(@(x) size(x.nodes,1),superagglos);
cumNumNodes = [0;cumsum(numNodes);sum(numNodes)];
edges = cell2mat(arrayfun(@(x) superagglos(x).edges+cumNumNodes(x),(1:numel(superagglos))','uni',0));% putting all edges together (inceasing index for each superagglo
edges = edges(all(ismember(edges,idsToKeep),2),:);

fNames = setdiff(fieldnames(superagglos),{'nodes','edges'});
tmpstrct = cell(numel(fNames),1);
for f = 1:numel(fNames)
    tmpstrct{f} = cat(1,superagglos.(fNames{f}));
end

if noCC
    aggloLUT = repelem((1:numel(superagglos))',numNodes)';
    newnodes = arrayfun(@(x) nodes(aggloLUT==x & toKeep,:) ,1:numel(superagglos),'uni',0);
    for f = 1:numel(fNames)
        if size(tmpstrct{f},1) == size(nodes,1) % apply same sorting to all other field names if their size was the same as the number of nodes, else leave the same
            tmpstrct{f} = arrayfun(@(x) tmpstrct{f}(aggloLUT==x & toKeep,:,:),1:numel(superagglos),'uni',0)';
        end
    end
    newedges = arrayfun(@(x) edges(all(edges > cumNumNodes(x) & edges <= cumNumNodes(x+1) ,2),:) ,1:numel(superagglos),'uni',0);
    emptyEdges = cellfun(@isempty,newedges);
    newedges(~emptyEdges) = cellfun(@(x) x - [sum(bsxfun(@ge,x(:,1),idsToDelete'),2),sum(bsxfun(@ge,x(:,2),idsToDelete'),2)],newedges(~emptyEdges),'uni',0);
    newSuperagglos = cell2struct([newedges;newnodes;tmpstrct(:)],[{'edges'},{'nodes'},fNames'],1);
    return
end




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
for f = 1:numel(fNames)
    if size(tmpstrct{f},1) == size(nodes,1) % apply same sorting to all other field names if their size was the same as the number of nodes, else leave the same
        tmpstrct{f} = cellfun(@(x) tmpstrct{f}(x,:),equivalenceClass,'uni',0)';
    end
end
% tranform the global node ids in the edge vector to local node ids
[~, newedges] = cellfun(@(x,y) ismember(x,y),newedges,equivalenceClass','uni',0);

newSuperagglos = cell2struct([newedges;newnodes;tmpstrct(:)],[{'edges'},{'nodes'},fNames'],1);
end

