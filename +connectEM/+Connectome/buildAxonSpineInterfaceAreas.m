function asiT = buildAxonSpineInterfaceAreas( ...
        param, graph, axons, shAgglos, synapses, synT)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    maxSegId = Seg.Global.getMaxSegId(param);
    axonLUT = Agglo.buildLUT(maxSegId, axons);
    shLUT = Agglo.buildLUT(maxSegId, shAgglos);
    
    synT.shId = cellfun( ...
        @(segIds) max(shLUT(segIds)), ...
        synapses.postsynId(synT.id));
    synT(~synT.shId, :) = [];
    
   [~, asiT, synIds] = unique(synT(:, {'preAggloId', 'shId'}), 'rows');
    asiT = synT(asiT, {'id', 'preAggloId', 'postAggloId', 'shId'});
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
