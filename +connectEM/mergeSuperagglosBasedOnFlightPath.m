function superagglos_new = mergeSuperagglosBasedOnFlightPath( ...
        superagglos, eqClassCCfull, startAgglo, endAgglo, ff)
    % Written by
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    %   Christian Schramm <christian.schramm@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % configuration
    voxelSize = [11.24, 11.24, 28];
    
    % sanity check
    % Make sure that flight paths actually do have a segment overlap with
    % their `startAgglo` and, if present, `endAgglo`.
    for curIdx = 1:numel(ff.filenames)
        curSegIdsHit = cat( ...
            2, ff.segIds{curIdx}, ff.neighbours{curIdx});
        curSegIdsHit = curSegIdsHit(curSegIdsHit > 0);
        curSegIdsHit = unique(curSegIdsHit);
        
        % make sure there is overlap with start agglomerate
        curStartSegIds = superagglos(startAgglo{curIdx}).nodes(:, 4);
        assert(numel(intersect(curSegIdsHit, curStartSegIds)) > 0);
        
        % same test if there is an end agglomerate
        if isempty(endAgglo{curIdx}); continue; end;
        curEndSegIds = superagglos(endAgglo{curIdx}).nodes(:, 4);
        assert(numel(intersect(curSegIdsHit, curEndSegIds)) > 0);
    end
    
    % find connected component for each flight path
    attachments = cellfun(@(x,y)[x y], startAgglo, endAgglo, 'uni', 0);
    
    ccLookup = Agglo.buildLUT(max(cell2mat(eqClassCCfull)), eqClassCCfull);
    attachmentsCC = cellfun(@(x)unique(ccLookup(x)), attachments);
    
    % Generate new superagglo
    for i=1:length(eqClassCCfull)
        % Concatenate superagglos of this equivalence class
        classLength = cell2mat(arrayfun(@(x)size(superagglos(x).nodes,1),eqClassCCfull{i},'uni',0));
%         classLookup = repelem([1:length(classLength)]',classLength);
        classLookup = repelem(eqClassCCfull{i},classLength);
        superagglos_new(i,1).nodes = cat(1, superagglos(eqClassCCfull{i}).nodes);
        
        % collect and renumber edges
        nrNodes = cumsum(arrayfun(@(x)size(x.nodes,1), superagglos(eqClassCCfull{i})));
        nodeOffset = cat(1, 0, nrNodes);
        edgeCell = arrayfun(@(x,y)x.edges+y, superagglos(eqClassCCfull{i}), nodeOffset(1:end-1), 'uni', 0);
        superagglos_new(i,1).edges = cat(1, edgeCell{:});
        
        % Find all queries that link to current superagglo (start or end)
        queryIdx = attachmentsCC == i;
        
        if any(queryIdx)
            % concatenate nodes from flight paths
            % by definition, flight paths nodes have `nan` as segment ID
            newNodes = cat(1, ff.nodes{queryIdx});
            newNodes(:,4) = NaN(size(newNodes,1),1);
            
            queryIndices = find(queryIdx);
            flightAgglos = arrayfun(@(x)cat(1,[startAgglo{x};endAgglo{x}]),queryIndices,'uni',0);
            % reconstruct flight path edges using MST
            edgeCell = cellfun(@(n) ...
                Graph.getMST(bsxfun(@times, n(:, 1:3), voxelSize)), ...
                ff.nodes(queryIdx), 'uni', 0);
            
            % collect and renumber edges
            nrNodes = cumsum(cellfun(@(x)size(x,1), ff.nodes(queryIdx)));
            nodeOffset = cat(2, 0, nrNodes);
            newEdges = cellfun(@(x,y)x+y, edgeCell, num2cell(nodeOffset(1:end-1)), 'uni', 0);
            newEdges = cat(1, newEdges{:});
            
            % Add in edges that connect current query to current superagglo
            % Just add one edge which is directly at the border between 
            % flight path and agglomerate
            segIdsThisSuperagglo = superagglos_new(i).nodes(:,4);
            segIdsHitByNewNodes = cat(2, cat(1, ff.segIds{queryIdx}), cat(1, ff.neighbours{queryIdx}));
            [idxNewNodes, idxOldNodes] = ismember(segIdsHitByNewNodes, segIdsThisSuperagglo);
            
            classOrigin = zeros(size(idxOldNodes));
            classOrigin(idxOldNodes>0) = classLookup(idxOldNodes(idxOldNodes>0));
            classOrigin = max(classOrigin')';
            classes = unique(classOrigin(classOrigin~=0));
            occurences = sum(idxNewNodes,2);
            

            cut = zeros(size(occurences))';
            N = numel(occurences);
            cut(nrNodes+1) = 1;
            cut(end) = [];
            ind = find(cut);
            ind_before = [ind-1 N]; ind_before(ind_before < 1) = 1;
            ind_after = [1 ind]; ind_after(ind_after > N) = N;
            occurences = arrayfun(@(x,y) occurences(x:y), ind_after, ind_before, 'uni', 0);
            classOrigin = arrayfun(@(x,y) classOrigin(x:y), ind_after, ind_before, 'uni', 0);

            occurences = cellfun(@(x,y)assignValues(x,y,1),occurences,occurences,'uni',0);
            classOrigin = cellfun(@(x,y)assignValues(x,y,1),classOrigin,occurences,'uni',0);
            
            % Eliminate evidences of not attached agglomerates
            eliminateClasses = cellfun(@(x,y)~ismember(x,y),classOrigin',flightAgglos,'uni',0);
            classOrigin = cellfun(@(x,y)assignValues(x,y),classOrigin,eliminateClasses','uni',0);

            nodesBeingAttached = cellfun(@(x)find(x),occurences,'uni',0);
            nodesClasses = cellfun(@(x,y)x(y),classOrigin,nodesBeingAttached,'uni',0);
            classes = cellfun(@(x)unique(x(x~=0)),classOrigin,'uni',0);

            borderNodeIdx = cell(size(occurences));
            for j=1:length(occurences)
                if numel(classes{j})==1
                    borderNode = max(find(nodesClasses{j}));
                    borderNode = nodesBeingAttached{j}(borderNode);
                    borderNodeIdx{j} = zeros(size(occurences{j}));
                    borderNodeIdx{j}(borderNode) = 1;
                elseif numel(classes{j})==2
                    firstNode = nodesClasses{j}(1);
                    currentNode=1;
                    while firstNode == nodesClasses{j}(currentNode)
                        currentNode = currentNode+1;
                    end
                    borderNodes = [currentNode-1 currentNode]';
                    borderNodes = nodesBeingAttached{j}(borderNodes);
                    borderNodeIdx{j} = zeros(size(occurences{j}));
                    borderNodeIdx{j}(borderNodes) = 1;
                elseif numel(classes{j})==0
                    disp('Warning: Not a single node has enough node evidence!')
                else
                    disp('Warning: One flight has enough node evidence to different end agglos!')
                end
            end
            
            borderNodeIdx = cat(1,borderNodeIdx{:});
            newNodesBeingAttached = find(borderNodeIdx);
            flatten = @(x)x(:);
            
            % Build segment-to-query edges
            connectingEdges = cat(2, flatten(idxOldNodes(newNodesBeingAttached,:)), repmat(newNodesBeingAttached + size(superagglos_new(i).nodes,1),27,1));
            connectingEdges = unique(connectingEdges, 'rows');
            connectingEdges(any(connectingEdges == 0,2),:) = [];
            
            % Add to output structure
            superagglos_new(i,1).edges = cat(1, superagglos_new(i).edges, newEdges + size(superagglos_new(i).nodes,1), connectingEdges);
            superagglos_new(i,1).nodes = cat(1, superagglos_new(i).nodes, newNodes);
        end
        
        % sanity checks
        % * edges are correctly sorted
        assert(all(diff(superagglos_new(i,1).edges, 1, 2) > 0));
        
        % * single connected component
        switch size(superagglos_new(i,1).nodes, 1)
            case 0
                assert(false);
            case 1
                assert(isempty(superagglos_new(i,1).edges));
            otherwise
                assert(numel(Graph.findConnectedComponents( ...
                    superagglos_new(i).edges)) == 1);
        end
    end
end

function input1 = assignValues(input1,input2,condition)
    if nargin < 3
        condition = false;
    end
    if condition
        input1(input2 < 14) = 0;
    else
        input1(input2) = 0;
    end
end