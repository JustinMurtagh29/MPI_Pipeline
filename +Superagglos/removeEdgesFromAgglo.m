function [newSuperagglos] = removeEdgesFromAgglo(superagglos,edgeIds,noCC )
% removes sets of edge indices from superagglos and performs connected
% component afterswards to ensure internode connectivity in each superagglo
if ~exist('noCC','var') || isempty(noCC)
    noCC = 0;
end

if ~ iscell(edgeIds)
    edgeIds = {edgeIds};
end
newSuperagglos = superagglos([]);
for n = 1:numel(superagglos)
%     nodes = cell2mat(arrayfun(@(x) x.nodes,superagglos(n),'uni',0));
    nodes = superagglos(n).nodes;
    fNames = setdiff(fieldnames(superagglos(n)),{'nodes','edges'});
    tmpstrct = cell(numel(fNames),1);
    for f = 1:numel(fNames)
        tmpstrct{f} = cat(1,superagglos(n).(fNames{f}));
    end   
    if islogical(edgeIds{n})
        toKeep = ~edgeIds{n};
    else
        toKeep = true(size(superagglos(n).edges,1),1);
        toKeep(edgeIds{n}) = false;
    end
  
    if noCC
        newSuperagglos(n) = superagglos(n);
        newSuperagglos(n).edges = superagglos(n).edges(toKeep,:);
        continue
    end
    edges = superagglos(n).edges(toKeep,:);
    % generate selfEdges of all superagglo segments
    selfEdges = repmat(unique(edges(:)),1,2);
    
    [equivalenceClass, aggloLUT] = Graph.findConnectedComponents(cat(1,selfEdges,edges),0,1);
    if isempty(equivalenceClass)
        continue
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
    newedges = cell(1,numel(equivalenceClass));
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
    [~, newedges] = cellfun(@(x,y) ismember(x,y),newedges,equivalenceClass','uni',0);
    
    newSuperagglos = cat(1,newSuperagglos,cell2struct([newedges;newnodes;tmpstrct{:}],[{'edges'},{'nodes'},fNames'],1));
end
end

