function asiT = buildAxonSpineInterfaces( ...
        param, graph, shAgglos, conn, syn, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opt = struct;
    opt.addBorderIdsVar = false;
    opt = Util.modifyStruct(opt, varargin{:});
    
    maxSegId = Seg.Global.getMaxSegId(param);
    axonLUT = Agglo.buildLUT(maxSegId, conn.axons);
    
    % NOTE(amotta): Find spine heads that overlap with axon agglomerates.
    % This is a sign that something went wrong. We don't want to use these
    % data points in our analysis.
    %   By removing the spine head we make sure that the corresponding
    % synapse is not detected as being onto a spine head, and thus won't
    % make it into the axon-spine interface table.
    shIds = find(cellfun(@(ids) ~any(axonLUT(ids)), shAgglos));
    shLUT = Agglo.buildLUT(maxSegId, shAgglos(shIds), shIds);
    
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
        graph.borderArea(mask), ...
       [height(asiT), 1], @sum, nan);
    
    if opt.addBorderIdsVar
        % NOTE(amotta): The `borderIds` variable is only addedn when
        % explicitly asked for by the user. It's rather slow to allocate so
        % many cells...
        asiT.borderIds = accumarray( ...
            graph.axonSpineId(mask), graph.borderIdx(mask), ...
           [height(asiT), 1], @(ids) {ids(:)}, {zeros(0, 1)});
    end
    
    % NOTE(amotta): The synapse agglomeration happened independently of the
    % neurite reconstruction. In very rare cases (55 times at the time of
    % writing) it can thus happen that there does not exist any border
    % between a axon and spine head agglomerates.
    %   For now, this case is eliminated by brute force.
    asiT(isnan(asiT.area), :) = [];
end
