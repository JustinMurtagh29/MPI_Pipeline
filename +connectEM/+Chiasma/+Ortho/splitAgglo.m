function newAgglos = splitAgglo(param, chiParam, agglo, queries)
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
    newEndings = cell(chiCount, 1);
    
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
            newEndings{curChiIdx} = zeros(0, 1);
            
            % mark chiasma as solved
            agglo.solvedChiasma(curNodeId) = true;
            continue;
        end
        
        %% mark new open endings
        dangExits = find(accumarray(curGrouping.groupId, 1) == 1);
        dangExits = curGrouping.id(ismember(curGrouping.groupId, dangExits));
        newEndings{curChiIdx} = curExitNodeIds(dangExits);

        %% cut out sphere
       [nodes, edges, nodeIds, ~, coreNodeIds] = ...
           connectEM.detectChiasmataPruneToSphere( ...
                param, agglo.nodesNm, agglo.edges, curNodeId);
        comps = Graph.findConnectedComponents( ...
            edges, false, false, numel(nodeIds));

        % drop tiny components
        isExit = cellfun(@(idx) max(pdist2( ...
            nodes(idx, :), curCenterNm)) > param.minNodeDist, comps);
        nonExitNodeIds = nodeIds(cell2mat(comps(~isExit)));
        
        % collect deleted nodes and edges
        nodesToDelete{curChiIdx} = cat( ...
            1, nonExitNodeIds(:), coreNodeIds);
        edgesToDelete{curChiIdx} = find(any( ...
            ismember(agglo.edges, nodesToDelete{curChiIdx}), 2));
        
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
    chiasmaValid = ~cellfun(@(e) any(ismember( ...
        e(:), cell2mat(nodesToDelete))), edgesToAdd);
    
    nodesToAdd = zeros(0, 4);
    nodesToDelete = cell2mat(nodesToDelete(chiasmaValid));
    edgesToDelete = cell2mat(edgesToDelete(chiasmaValid));
    edgesToAdd = cell2mat(edgesToAdd(chiasmaValid));
    newEndings = cell2mat(newEndings(chiasmaValid));
    
    % build mask for endings
    agglo.endingMask = false(size(agglo.nodes, 1), 1);
    agglo.endingMask(agglo.endings) = true;
    agglo.endingMask(newEndings) = true;
    agglo = rmfield(agglo, 'endings');
    
    patchExtra = struct;
    patchExtra.endingMask = false(size(nodesToAdd, 1), 1);
    patchExtra.solvedChiasma = false(size(nodesToAdd, 1), 1);
    
    % split axon
    newAgglos = Superagglos.patch( ...
        agglo, nodesToAdd, nodesToDelete, ...
        edgesToAdd, edgesToDelete, patchExtra);
    
	% build endings
    newAggloEndings = arrayfun( ...
        @(a) find(a.endingMask(:)), ...
        newAgglos, 'UniformOutput', false);
   [newAgglos.endings] = deal(newAggloEndings{:});
    clear newAggloEndings;
    
    newAgglos = rmfield(newAgglos, 'endingMask');
    newAgglos = orderfields(newAgglos, { ...
        'nodes', 'edges', 'endings', 'solvedChiasma'});
end