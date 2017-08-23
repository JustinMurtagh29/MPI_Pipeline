function [ newSuperagglos, equivalenceClass ] = executeEdges(superagglos, edgesToExecute, segmentMeta)
% executes the given edges between superagglos, and also adds
% single segments of these edges

% Extract all segment ID based edges from superagglos
aggloEdges = cell2mat(arrayfun(@(x) reshape(x.nodes(x.edges,4),[],2),...
    superagglos,'uni',0));
if all(ismember(edgesToExecute,aggloEdges,'rows'))
    % In case all edgesToExecute are already present in superagglos
    newSuperagglos = superagglos;
    equivalenceClass = num2cell(1:numel(superagglos));
else
    % get all seg ids that are part of superagglos including single segment
    % superagglos
    segIds = cell2mat(arrayfun(@(x)x.nodes(:,4), superagglos, 'uni', 0));
    
    % Concatenate all agglo edges and edges to be executed
    allEdges = cat(1,aggloEdges, edgesToExecute);
    
    % Connected components on all edges
    [equivalenceClass, aggloLUT] = Graph.findConnectedComponents(allEdges,0,0);
    
    % create boolean which equivalence classes contain single edges
    singleSegAgglos = true(numel(equivalenceClass),1);
    singleSegAgglos(aggloLUT(allEdges(:))) = false;
    
    % get indices to equivalence classes which have segments that were part
    % of the superagglo. These are kept later.
    eqClassesToKeep = unique(aggloLUT(segIds));
    
    % sort hot edges according to the agglos they will be put in
    [saggloLUT,sIdx] = sort(aggloLUT);
    allEdges = allEdges(sIdx,:);
    
    % transform the edges into cell
    newedges = cell(numel(aggloOld),1);
    % Treat agglomerates containing single segment separately
    newedges(~singleSegAgglos) = mat2cell(allEdges,histc(saggloLUT,unique(saggloLUT)));
    newedges(singleSegAgglos) = {zeros(0,2)}; % give empty edges correct dimension
    
    % filter edge and equivalence class list to keep only those which have
    % parts of the original superagglo
    newedges = newedges(eqClassesToKeep);
    equivalenceClass = equivalenceClass(eqClassesToKeep);
    
    %  create node cell array including segID information
    newnodes = cellfun(@(x) [segmentMeta.point(:,x)',x],equivalenceClass);

    newSuperagglos = cell2struct([newedges';newnodes'],{'edges','nodes'},1);
end

% Check whether classes were mutually exclusive
assert(numel(segIds) == numel(unique(segIds)));

end

