function newAgglos = splitWithQueries(param, chiParam, agglo, queries)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % update parameter structure
    kvPairs = cat(2, fieldnames(chiParam), struct2cell(chiParam));
    kvPairs = transpose(kvPairs);
    
    param = Util.modifyStruct(param, kvPairs{:});
    clear kvPairs;
    
    %% split chiasmata
    agglo.nodesNm = agglo.nodes(:, 1:3) .* param.raw.voxelSize;
    
    % generate fake edge probabilities
    % These are not actually of any importance.
    edgeProb = ones(size(agglo.edges, 1), 1);
    nodeCount = size(agglo.nodes, 1);
    
    % state
    chiCount = size(queries, 1);
    nodesToDelete = cell(chiCount, 1);
    edgesToDelete = cell(chiCount, 1);
    
    edgesToAdd = cell(chiCount, 1);
    edgeBuilder = @(nodeIds) cat(2, ...
        repelem(nodeIds(1), numel(nodeIds) - 1, 1), ...
        reshape(nodeIds(2:end), [], 1));
    
    for curChiIdx = 1:chiCount
        curChi = queries(curChiIdx, :);
        
        curNodeId = curChi.centerNodeId;
        curCenterNm = agglo.nodesNm(curNodeId, :);

        %% cut out sphere
       [nodes, edges, ~, nodeIds] = ...
           connectEM.detectChiasmataPruneToSphere( ...
                agglo.nodesNm, agglo.edges, edgeProb, param, curNodeId);    
        comps = Graph.findConnectedComponents(edges, false);

        % drop tiny components
        isExit = cellfun(@(idx) max(pdist2( ...
            nodes(idx, :), curCenterNm)) > param.minNodeDist, comps);
      
        nrExits = sum(isExit);
        nonExitNodeIds = nodeIds(cell2mat(comps(~isExit)));
        
        % collect deleted nodes and edges
        nodesToDelete{curChiIdx} = cat( ...
            1, nonExitNodeIds(:), setdiff( ...
            reshape(1:nodeCount, [], 1), nodeIds));
        edgesToDelete{curChiIdx} = find( ...
            any(ismember(agglo.edges, nonExitNodeIds), 2) ...
         | ~ismember(agglo.edges, nodeIds(edges), 'rows'));
        
        %% build new edges
        curGrouping = curChi.exits{1};
        curExitNodeIds = curChi.exitNodeIds{1};
     
        % sanity check
        if ~(isequal(curGrouping.id, transpose(1:nrExits)))
            keyboard
        end
        
        % build edges to add
        curEdgesToAdd = accumarray( ...
            curGrouping.groupId, curExitNodeIds, [], ...
            @(nodeIds) {edgeBuilder(nodeIds)}, {zeros(0, 2)});
        curEdgesToAdd = cat(1, curEdgesToAdd{:});
        curEdgesToAdd = sort(curEdgesToAdd, 2);
        
        edgesToAdd{curChiIdx} = curEdgesToAdd;
    end
    
    %% patch in new edges
    chiasmaValid = ~cellfun(@(e) any( ...
        ismember(e(:), cell2mat(nodesToDelete))), edgesToAdd);
    
    nodesToDelete = cell2mat(nodesToDelete(chiasmaValid));
    edgesToDelete = cell2mat(edgesToDelete(chiasmaValid));
    edgesToAdd = cell2mat(edgesToAdd(chiasmaValid));
    
    % TODO(amotta):
    % Most of the code below was copied in from splitChiasmataMulti and
    % should be factored out into a reusable function!
    nodesToKeep = setdiff(1:size(agglo.nodes, 1), nodesToDelete);
    newNodes = agglo.nodes(nodesToKeep, :);
    
    newEdges = cat(1, agglo.edges, edgesToAdd);
    newEdges(edgesToDelete, :) = [];
   [~, newEdges] = ismember(newEdges, nodesToKeep);
    assert(all(newEdges(:)));
   
    newNodeComps = Graph.findConnectedComponents(newEdges);
    newNodeComps = cat(1, newNodeComps, num2cell(reshape( ...
        setdiff(size(newNodes, 1), cell2mat(newNodeComps)), [], 1)));
    
    newNodeLUT = Agglo.buildLUT(size(newNodes, 1), newNodeComps);
    newEdgeComps = newNodeLUT(newEdges(:, 1));
    
    newAggloCount = numel(newNodeComps);
    newAgglos = struct;
    
    for curIdx = 1:newAggloCount
        curNodeIds = newNodeComps{curIdx};
        curNodes = newNodes(curNodeIds, :);
        % curSolvedChiasma = newSolvedChiasma(curNodeIds);
        % curEndings = find(ismember(curNodeIds, newEndings));
        
        curEdges = newEdges(newEdgeComps == curIdx, :);
        [~, curEdges] = ismember(curEdges, curNodeIds);
        assert(all(curEdges(:)));
        
        newAgglos(curIdx).nodes = curNodes;
        newAgglos(curIdx).edges = curEdges;
        % newAgglos(curIdx).endings = curEndings;
        % newAgglos(curIdx).solvedChiasma = curSolvedChiasma;
    end
    
    newAgglos = reshape(newAgglos, [], 1);
end