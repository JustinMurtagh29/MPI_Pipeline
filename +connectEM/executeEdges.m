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
%     % Generate self edges for each segId (necessary for
%     % findConnectedComponents on single segment agglos)
%     segIds = cell2mat(arrayfun(@(x)x.nodes(:,4), superagglos, 'uni', 0));
%     selfEdges = repmat(segIds, 1, 2);
    segIds = cell2mat(arrayfun(@(x)x.nodes(:,4), superagglos, 'uni', 0));
    % Concatenate all self, agglo and edges to be executed
    allEdges = cat(1,aggloEdges, edgesToExecute);
    % Connected components on all edges
    [equivalenceClass, aggloLUT] = Graph.findConnectedComponents(allEdges,0,0);
    
    singleSegAgglos = true(numel(equivalenceClass),1);
    singleSegAgglos(allEdges(:)) = false;
    
    eqClassesToKeep = unique(aggloLUT(segIds));
%     edgesAggloBased = aggloLUT(allEdges);
    
    % sort hot edges according to the agglos they will be put in
    [saggloLUT,sIdx] = sort(aggloLUT);
    allEdges = allEdges(sIdx,:);
    
    % transform the edges into cell and create node cell array including segID
    % information
    newedges = cell(numel(aggloOld),1);
    % Treat agglomerates containing single segment separately
    newedges(~singleSegAgglos) = mat2cell(allEdges,histc(saggloLUT,unique(saggloLUT)));
    newedges(singleSegAgglos) = {zeros(0,2)}; % give empty edges correct dimension
    
    newedges = newedges(eqClassesToKeep);
    equivalenceClass = equivalenceClass(eqClassesToKeep);
    
    %
    %     [~,idx] = setdiff(1:numel(aggloLUT),segIds);
%     aggloLUT(idx) = [];
%     idxToKeep = unique(aggloLUT);
%     idxToDelete = true(numel(equivalenceClass),1);

    %     idxToDelete(idxToKeep) = false;
%     % Now remove everything that does not contain an inital superagglo
%     equivalenceClass = equivalenceClass(idxToKeep);
%     
%     aggloLUT = connectEM.changeem(aggloLUT, (1:numel(equivalenceClass)) - cumsum(idxToDelete),1:numel(equivalenceClass))
    
    newnodes = cellfun(@(x) [segmentMeta.point(:,x)',x],equivalenceClass);
%     [~, newedges] = cellfun(@(x) ismember(allEdges(aggloLUT(allEdges(:,1))==x,:),newnodes{x}(:,4)),1:numel(equivalenceClass),'uni',0);
    newSuperagglos = cell2struct([newedges';newnodes'],{'edges','nodes'},1);
    
end

% Check whether classes were mutually exclusive
assert(numel(segIds) == numel(unique(segIds)));

end

