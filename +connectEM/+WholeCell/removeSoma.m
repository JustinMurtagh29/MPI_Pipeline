function agglo = removeSoma(param, agglo, soma)
    % agglo = removeSoma(param, agglo, soma)
    %   Removes all soma-overlapping segments from a super-agglomerate. If
    %   the super-agglomerate is split into multiple connected components,
    %   the components are joined via a non-segment node placed at the
    %   center-of-mass of the soma.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % Find node IDs to keep
    nodeIds = Agglo.fromSuperAgglo(soma);
    nodeIds = find(~ismember(agglo.nodes(:, 4), nodeIds));
    
    % Drop soma nodes
    agglo.nodes = agglo.nodes(nodeIds, :);
   [~, agglo.edges] = ismember(agglo.edges, nodeIds);
    agglo.edges(~all(agglo.edges, 2), :) = [];
    
    % Clean up edges
    agglo = SuperAgglo.clean(agglo, false);
    nodeCount = size(agglo.nodes, 1);
    
    % Check if super-agglomerate is still connected
    graph = sparse( ...
        nodeCount, nodeCount, true, ...
        agglo.edges(:, 2), agglo.edges(:, 1));
    graph(1, 1) = true;
    
   [compCount, compIds] = graphconncomp(graph, 'Directed', false);
    if compCount < 2; return; end
    
    % Nope, there are multiple components. Let's fix that
    warning('Connecting components via virtual soma node');
    
    voxelSize = param.raw.voxelSize;
    somaPos = mean(soma.nodes(:, 1:3), 1);
    somaPosNm = somaPos .* voxelSize;
    
    % Find node closest to soma for each component
    compNodeIds = accumarray( ...
        compIds, (1:nodeCount)', ...
        [], @(nodeIds) {nodeIds});
   	compNodeIds = cellfun( ...
        @(ids) ids(secondOutput(@() pdist2( ...
            voxelSize .* agglo.nodes(ids, 1:3), ...
            somaPosNm, 'euclidean', 'Smallest', 1))), ...
        compNodeIds);
    
    % Add new node and edges
    newNode = [round(somaPos), nan];
    newEdges = nan(compCount, 2);
    newEdges(:, 1) = compNodeIds;
    newEdges(:, 2) = nodeCount + 1;
    
	agglo.nodes = [agglo.nodes; newNode];
    agglo.edges = [agglo.edges; newEdges];
    agglo = SuperAgglo.clean(agglo);
end

function b = secondOutput(f)
    [~, b] = f();
end
