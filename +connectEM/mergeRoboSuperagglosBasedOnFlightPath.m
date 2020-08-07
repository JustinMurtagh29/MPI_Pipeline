function superagglos_new = mergeRoboSuperagglosBasedOnFlightPath( ...
        superagglos, eqClassCCfull, startAgglo, endAgglo, ff,allowMultiAttachments, checkOccurences, useNeighbours, assertDoublyAttached)
    % INPUTS 
    % allowMultiAttachments: Boolean allowing attachment of the flight paths with
    %                        multiple ends (DEFAULT 0)
    % checkOccurences: Boolean to (re)check node edvidence of flight paths
    %                  with segments (DEFAULT 1). Should be switched to
    %                  False if flight tracings were done on a different
    %                  agglomerate state than the one they are now attached
    %                  to to avoid problems (node evidence was checked
    %                  anyway before).
    % useNeighbours:   Whether to use ff.neighbours information about segmentIds hit
    %                  within the 26 neighborhood of a voxel. Usually, we don't need
    %                  to do that for RoboEM.
    % assertDoublyAttached: activates a few assert statements that should hold, when
    %                       each tracing should always attach to two agglos
    %
    % Originally written by
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    %   Christian Schramm <christian.schramm@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    %   Marcel Beining <marcel.beining@brain.mpg.de>
    % Copied from connectEM.mergeSuperagglosBasedOnFlightPath and adapted to RoboEM usecase by
    %   Martin Schmidt <martin.schmidt@brain.mpg.de>
    
    % configuration
    voxelSize = [11.24, 11.24, 28];
    if ~exist('allowMultiAttachments','var') || isempty(allowMultiAttachments)
        allowMultiAttachments = false;
    end
    if ~exist('checkOccurences','var') || isempty(checkOccurences)
        checkOccurences = false;
    end
    if ~exist('useNeighbours', 'var') || isempty(useNeighbours)
        useNeighbours = false;
    end
    if ~exist('assertDoublyAttached', 'var') || isempty(assertDoublyAttached)
        assertDoublyAttached = true;
    end
    
    % sanity check
    % Make sure that flight paths actually do have a segment overlap with
    % their `startAgglo` and, if present, `endAgglo`.
    for curIdx = 1:numel(ff.segIds)
        if useNeighbours
            curSegIdsHit = cat( ...
                2, ff.segIds{curIdx}, ff.neighbours{curIdx});
        else
            curSegIdsHit = ff.segIds{curIdx};
        end
        curSegIdsHit = curSegIdsHit(curSegIdsHit > 0);
        curSegIdsHit = unique(curSegIdsHit);
        
        % make sure there is overlap with start agglomerate
        curStartSegIds = superagglos(startAgglo{curIdx}).nodes(:, 4);
        assert(~isempty(intersect(curSegIdsHit, curStartSegIds)));
        
        % same test if there is an end agglomerate
        if isempty(endAgglo{curIdx}); continue; end
        for curEndIdx = 1:numel(endAgglo{curIdx})
            curEndSegIds = superagglos(endAgglo{curIdx}(curEndIdx)).nodes(:, 4);
            assert(~isempty(intersect(curSegIdsHit, curEndSegIds)));
        end
    end
    
    % find connected component for each flight path
    attachments = cellfun(@(x,y)[x; y], startAgglo, endAgglo, 'uni', 0);
    
    ccLookup = Agglo.buildLUT(max(cell2mat(eqClassCCfull)), eqClassCCfull); % creates lookup which agglo is in which eqClass
    attachmentsCC = cellfun(@(x)unique(ccLookup(x)), attachments);  % tells which flight path belongs to which eqClass. if a flight path would belong to several eqClasses, this would give an error
    
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
            
            flightAgglos = attachments(queryIdx);  % agglos that should be connected with the flight path(s)
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
            if useNeighbours
                segIdsHitByFlightPath = cat(2, cat(1, ff.segIds{queryIdx}), cat(1, ff.neighbours{queryIdx}));
            else
                segIdsHitByFlightPath = cat(1, ff.segIds{queryIdx});
            end
            [idxNewNodes, idxOldNodes] = ismember(segIdsHitByFlightPath, segIdsThisSuperagglo);
            
            classOrigin = zeros(size(idxOldNodes));
            classOrigin(idxOldNodes>0) = classLookup(idxOldNodes(idxOldNodes>0));
            % classOrigin = max(classOrigin,[],2); % > this does not make ANY sense and screws up things!! Why take the larger aggloId per node instead of taking the one that has higher evidence?!
            classOriginAsCell = num2cell(classOrigin, 2);
            classOriginPerNodeUq = cellfun(@unique, classOriginAsCell, 'uni', false);
            classOriginMaxCountPerNodeUq = cellfun(@(cO, cOUq) cOUq(argmaxfirst(sum(cOUq(:)' == cO(:), 2))), classOriginAsCell, classOriginPerNodeUq, 'uni', false);
            classOriginEliminate = cellfun(@(cO, cOMax) cO ~= cOMax, classOriginAsCell, classOriginMaxCountPerNodeUq, 'uni', false);
            classOriginEliminate = vertcat(classOriginEliminate{:});

            classOrigin = vertcat(classOriginMaxCountPerNodeUq{:});
            idxNewNodes(classOriginEliminate) = 0;
            idxOldNodes(classOriginEliminate) = 0;

            occurences = sum(idxNewNodes,2);
            

            % cut = zeros(size(occurences))';
            N = numel(occurences);
            % cut(nrNodes+1) = 1;
            % cut(end) = [];
            % ind = find(cut);            

            % assert(all(ind == nrNodes(1:end-1)+1)); % (note mschmidt: if this holds true, then creating cut is kind of unnecessary?!)

            ind = nrNodes(1:end-1)+1;

            ind_before = [ind-1 N]; ind_before(ind_before < 1) = 1;
            ind_after = [1 ind]; ind_after(ind_after > N) = N;
            occurences = arrayfun(@(x,y) occurences(x:y), ind_after, ind_before, 'uni', 0);
            classOrigin = arrayfun(@(x,y) classOrigin(x:y), ind_after, ind_before, 'uni', 0);
            
            if checkOccurences
                occurences = cellfun(@(x,y)assignValues(x,y,1),occurences,occurences,'uni',0);
                classOrigin = cellfun(@(x,y)assignValues(x,y,1),classOrigin,occurences,'uni',0);
            end
            % Eliminate evidences of not attached agglomerates
            eliminateClasses = cellfun(@(x,y)~ismember(x,y),classOrigin',flightAgglos,'uni',0);
            occurences = cellfun(@(x,y)assignValues(x,y),occurences,eliminateClasses','uni',0);
            classOrigin = cellfun(@(x,y)assignValues(x,y),classOrigin,eliminateClasses','uni',0);

            nodesBeingAttached = cellfun(@(x)find(x),occurences,'uni',0); % index to node ids of FLIGHT PATH (=newNodes) that have nodeEvidence
            if assertDoublyAttached
                assert(all(numel(nodesBeingAttached) > 0));
            end
            nodesClasses = cellfun(@(x,y)x(y),classOrigin,nodesBeingAttached,'uni',0);
            classes = cellfun(@(x)unique(x(x~=0)),classOrigin,'uni',0);

            borderNodeIdx = cell(size(occurences));
            for j=1:length(occurences)
                if assertDoublyAttached
                    assert(numel(classes{j}) == 2);
                end
                switch numel(classes{j})
                    case 1
                        borderNodes = find(nodesClasses{j},1,'last');
                        borderNodes = nodesBeingAttached{j}(borderNodes);
                    case 2
                        firstNode = nodesClasses{j}(1);
                        currentNode=1;
                        while firstNode == nodesClasses{j}(currentNode)
                            currentNode = currentNode+1;
                        end
                        borderNodes = [currentNode-1 currentNode]';
                        borderNodes = nodesBeingAttached{j}(borderNodes);
                    case 0
                        warning('Not a single node has enough node evidence!')
                        borderNodes = [];
                    otherwise
                        if allowMultiAttachments
                            borderNodes = find(abs(diff(nodesClasses{j},1,1)) > 0)+1;
                            borderNodes = cat(1,borderNodes-1,borderNodes);
                        else
                            warning('One flight has enough node evidence to different end agglos!')
                            borderNodes = [];
                        end
                end
                borderNodeIdx{j} = zeros(size(occurences{j}));
                borderNodeIdx{j}(borderNodes) = 1;
            end
            newNodesBeingAttached = find(cat(1,borderNodeIdx{:}));
            if assertDoublyAttached
                assert(length(newNodesBeingAttached) == 2*sum(queryIdx));
            end
            flatten = @(x)x(:);
            
            % Build segment-to-query edges
            connectingEdges = cat(2, flatten(idxOldNodes(newNodesBeingAttached,:)), repmat(newNodesBeingAttached + size(superagglos_new(i).nodes,1),1+26*useNeighbours,1));
            connectingEdges = unique(connectingEdges, 'rows');
            connectingEdges(any(connectingEdges == 0,2),:) = [];
            if assertDoublyAttached
                assert(size(connectingEdges, 1) >= 2*sum(queryIdx));
            end
            
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

function idx = argmaxfirst(a)
   [~, idx] = max(a);
   idx = idx(1);
end
