function mask = detectThalamocorticals(conn, interSyn)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    %% calculate inter-synapse distances
    % Previously we've calculate the synapse-to-synapse distances along the
    % axon. To calculate the inter-synapse distances we now need to construct a
    % tree-like representation of the axon. Let's do this by calculating the
    % minimal spanning-tree over the synapses.
    conn.axonMeta.pathLen = nan(size(conn.axonMeta.id));
    conn.axonMeta.interSynDists = cell(size(conn.axonMeta.id));

    for curIdx = 1:numel(interSyn.axonIds)
        curAxonId = interSyn.axonIds(curIdx);
        curPathLen = interSyn.axonPathLens(curIdx) / 1E3;

        curSynToSynDists = interSyn.synToSynDists{curIdx};
        curSynToSynDists = curSynToSynDists ./ 1E3;

        % NOTE(amotta): Zeros in the adjacency matrix are interpreted as
        % missing edge. This is a problem since synapse-to-synapse distance
        % zero occurs naturally when multiple synapses originate from the same
        % segment. Let's instead set these zeros to the smallest possible
        % non-zero value.
        curSynToSynDists(~curSynToSynDists) = eps;

        % claculate inter-synapse distances
        curInterSynDists = graphminspantree(sparse(curSynToSynDists));
        curInterSynDists = nonzeros(curInterSynDists);

        conn.axonMeta.pathLen(curAxonId) = curPathLen;
        conn.axonMeta.interSynDists{curAxonId} = curInterSynDists;
    end

    %% orthogonal approach
    minSynCount = 10;

    % empirically found
    minSynDensity = 0.25;
    minSpineSynFrac = 0.5;
    maxMedianInterSynDist = 2.2;

    % spine synapse fraction
    conn.axonMeta.spineFrac = ...
        conn.axonMeta.spineSynCount ...
        ./ conn.axonMeta.synCount;

    % synapse density
    conn.axonMeta.synDensity = ...
        conn.axonMeta.synCount ...
        ./ conn.axonMeta.pathLen;

    % median inter-synapse distance
    conn.axonMeta.medianInterSynDist = ...
        cellfun(@median, conn.axonMeta.interSynDists);

    mask = ...
        (conn.axonMeta.synCount >= minSynCount) ...
      & (conn.axonMeta.spineFrac >= minSpineSynFrac) ...
      & (conn.axonMeta.synDensity >= minSynDensity) ...
      & (conn.axonMeta.medianInterSynDist <= maxMedianInterSynDist);
end