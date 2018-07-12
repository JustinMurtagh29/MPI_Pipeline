function synIds = getSynapses(conn)
    % synIds = getSynapses(conn)
    %   Returns a cell array, where `synIds{i}` contains the indices of the
    %   synapses of axon `i`.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    connT = conn.connectome;
   [~, synIds] = ismember( ...
        connT.edges(:, 1), conn.axonMeta.id);
    synIds = accumarray( ...
        synIds, (1:size(synIds, 1))', size(conn.axons), ...
        @(ids) {sort(cell2mat(connT.synIdx(ids)))}, {zeros(0, 1)});
    
    % Sanity check
    assert(isequal(cellfun(@numel, synIds), conn.axonMeta.synCount));
end
