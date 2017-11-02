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
        
        %% check if anything to do
        curGrouping = curChi.exits{1};
        curExitNodeIds = curChi.exitNodeIds{1};
        
        % check if fully connected
        if max(curGrouping.groupId) == 1
            % no modifications needed
            nodesToDelete{curChiIdx} = zeros(0, 1);
            edgesToDelete{curChiIdx} = zeros(0, 1);
            edgesToAdd{curChiIdx} = zeros(0, 2);
            
            % mark chiasma as solved
            agglo.solvedChiasma(curNodeId) = true;
            continue;
        end

        %% cut out sphere
       [nodes, edges, ~, nodeIds] = ...
           connectEM.detectChiasmataPruneToSphere( ...
                agglo.nodesNm, agglo.edges, edgeProb, param, curNodeId);    
        comps = Graph.findConnectedComponents(edges, false);

        % drop tiny components
        isExit = cellfun(@(idx) max(pdist2( ...
            nodes(idx, :), curCenterNm)) > param.minNodeDist, comps);
        nonExitNodeIds = nodeIds(cell2mat(comps(~isExit)));
        
        % collect deleted nodes and edges
        nodesToDelete{curChiIdx} = cat( ...
            1, nonExitNodeIds(:), setdiff( ...
            reshape(1:nodeCount, [], 1), nodeIds));
        edgesToDelete{curChiIdx} = find( ...
            any(ismember(agglo.edges, nonExitNodeIds), 2) ...
         | ~ismember(agglo.edges, nodeIds(edges), 'rows'));
        
        %% build new edges
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
    
    nodesToAdd = zeros(0, 4);
    nodesToDelete = cell2mat(nodesToDelete(chiasmaValid));
    edgesToDelete = cell2mat(edgesToDelete(chiasmaValid));
    edgesToAdd = cell2mat(edgesToAdd(chiasmaValid));
    
    newAgglos = Superagglos.patch( ...
        agglo, nodesToAdd, nodesToDelete, edgesToAdd, edgesToDelete);
end