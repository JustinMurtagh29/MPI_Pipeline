function asiT = buildAxonSpineInterfaces( ...
        param, graph, shAgglos, conn, syn)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    maxSegId = Seg.Global.getMaxSegId(param);
    axonLUT = Agglo.buildLUT(maxSegId, conn.axons);
    shLUT = Agglo.buildLUT(maxSegId, shAgglos);
    
    synT = connectEM.Connectome.buildSynapseTable(conn, syn);
    synT.type = syn.synapses.type(synT.id);
    
    synT.shId = cellfun( ...
        @(segIds) max(shLUT(segIds)), ...
        syn.synapses.postsynId(synT.id));
    synT(~synT.shId, :) = [];
    
   [~, asiT, synIds] = unique(synT(:, {'preAggloId', 'shId'}), 'rows');
    asiT = synT(asiT, {'id', 'type', 'preAggloId', 'postAggloId', 'shId'});
    asiT.synIds = accumarray(synIds, synT.id, [], @(ids) {ids(:)});
    
    graph.shId = max(shLUT(graph.edges), [], 2);
    graph.preAggloId = max(axonLUT(graph.edges), [], 2);
    
   [~, graph.axonSpineId] = ismember( ...
        graph(:, {'preAggloId', 'shId'}), ...
        asiT(:, {'preAggloId', 'shId'}), 'rows');
    mask = (graph.axonSpineId ~= 0);
    
    asiT.area = accumarray( ...
        graph.axonSpineId(mask), ...
        graph.borderArea(mask));
end
