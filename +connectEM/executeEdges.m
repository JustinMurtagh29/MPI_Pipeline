function [ newSuperagglos, equivalenceClass ] = executeEdges(superagglos, edgesToExecute, segmentMeta)
% executes the given edges between superagglos, and also adds
% single segments of these edges

% get all seg ids that are part of superagglos including single segment
% superagglos
segIds = cell2mat(arrayfun(@(x)x.nodes(:,4), superagglos, 'uni', 0));
    
% check for duplicates in the segIds vector. This means there were overlaps
% between superagglos (or at least intra-agglo duplicates)
if numel(segIds) ~= numel(unique(segIds))
    warning('Detected overlapping segments in agglos. These will be merged, too, when re-executing all edges now!')
end

% Extract all segment ID based edges from superagglos
aggloEdges = cell2mat(arrayfun(@(x) reshape(x.nodes(x.edges,4),[],2),...
    superagglos,'uni',0));
if all(ismember(edgesToExecute,aggloEdges,'rows'))
    % In case all edgesToExecute are already present in superagglos
    newSuperagglos = superagglos;
    equivalenceClass = num2cell(1:numel(superagglos));
else
    % generate selfEdges of all superagglo segments
    selfEdges = repmat(segIds,1,2);
    
    % Concatenate all agglo edges and edges to be executed and remove
    % duplicates
    allNonSelfEdges = cat(1,aggloEdges, edgesToExecute);
    allNonSelfEdges = unique(sort(allNonSelfEdges,2),'rows');

    % Connected components on all edges including self edges
    [equivalenceClass, aggloLUT] = Graph.findConnectedComponents(cat(1,allNonSelfEdges,selfEdges),0,0);

    % create boolean which equivalence classes contain single edges
    singleSegAgglos = true(numel(equivalenceClass),1);
    singleSegAgglos(aggloLUT(allNonSelfEdges(:))) = false;
    
    % get indices to equivalence classes which have segments that were part
    % of the superagglo. These are kept later.
    eqClassesToKeep = unique(aggloLUT(segIds));
    
    % create lookup table which edge of allNonSelfEdges belongs to which
    % equivalence class
    edgesLUT = aggloLUT(allNonSelfEdges(:,1));
    
    % sort this LUT according to the agglos they will be put in and sort
    % edges accordingly
    [saggloLUT,sIdx] = sort(edgesLUT);
    allNonSelfEdges = allNonSelfEdges(sIdx,:);

    % transform the edges into cell
    newedges = cell(numel(superagglos),1);
    % Treat agglomerates containing single segment separately
    newedges(~singleSegAgglos) = mat2cell(allNonSelfEdges,histc(saggloLUT,unique(saggloLUT)));
    newedges(singleSegAgglos) = {zeros(0,2)}; % give empty edges correct dimension
       
    % filter edge and equivalence class list to keep only those which have
    % parts of the original superagglo
    newedges = newedges(eqClassesToKeep)';
    equivalenceClass = equivalenceClass(eqClassesToKeep);
    
    %  create node cell array including segID information
    newnodes = cellfun(@(x) [segmentMeta.point(:,x)',x],equivalenceClass,'uni',0)';

    % tranform the segIds in the edge vector to index to the node
    [~, newedges] = cellfun(@(x,y) ismember(x,y(:,4)),newedges,newnodes,'uni',0);

    newSuperagglos = cell2struct([newedges;newnodes],{'edges','nodes'},1);
end

end

