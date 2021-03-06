function [synIds, inConn] = getSynapses(conn, syn)
    % [synIds, inConn] = getSynapses(conn, syn)
    %   Returns a cell array, where `synIds{i}` contains the indices of the
    %   synapses of axon `i`. `inConn{i}(j)` is a logical indicating
    %   whether synapse `synIds{i}(j)` is also contained in the connectome.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    maxSegId = max( ...
        double(max(cellfun(@max, conn.axons))), ...
        double(max(cellfun(@max, syn.synapses.presynId))));
    axonLUT = Agglo.buildLUT(maxSegId, conn.axons);
    
    inConn = false(height(syn.synapses), 1);
    inConn(cell2mat(conn.connectome.synIdx)) = true;
    
    axonIds = cellfun( ...
        @(ids) reshape(setdiff(axonLUT(ids), 0), [], 1), ...
        syn.synapses.presynId, 'UniformOutput', false);
    
    synIds = repelem( ...
        reshape(1:height(syn.synapses), [], 1), ...
        cellfun(@numel, axonIds));
    axonIds = cell2mat(axonIds);
    
    synIds = accumarray( ...
        axonIds, synIds, size(conn.axons), ...
        @(ids) {reshape(ids, [], 1)}, {zeros(0, 1)});
    inConn = cellfun( ...
        @(ids) inConn(ids), synIds, ...
        'UniformOutput', false);
end
