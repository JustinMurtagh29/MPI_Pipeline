function [axonClasses, conn] = buildAxonClasses(conn, varargin)
    % axonClasses = buildAxonClasses(conn, varargin)
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opts = struct;
    opts.minSynPre = 1;
    opts.minSpineSynDensTc = 0.19;
    opts = Util.modifyStruct(opts, varargin{:});
    
    conn.axonMeta.fullPriSpineSynFrac = ...
        conn.axonMeta.fullPriSpineSynCount ...
     ./ conn.axonMeta.fullSynCount;
    conn.axonMeta.fullPriSpineSynDens = ...
        conn.axonMeta.fullPriSpineSynCount ...
     ./ conn.axonMeta.pathLen;
    
    % Local copy with modifications
    axonMeta = conn.axonMeta;
    axonMeta.fullSecSpineSynFrac = ...
        axonMeta.fullSpineSynCount ...
      - axonMeta.fullPriSpineSynCount;
    axonMeta.fullSecSpineSynFrac = ...
        axonMeta.fullSecSpineSynFrac ...
     ./ axonMeta.fullSynCount;
    
    excSomaMask = find( ...
       ~conn.denMeta.isInterneuron ...
      & conn.denMeta.targetClass == 'Somata');
    excSomaMask = ismember( ...
        conn.connectome.edges(:, 2), excSomaMask);
    
    axonMeta.excSomaSynCount = accumarray( ...
        conn.connectome.edges(excSomaMask, 1), ...
        cellfun(@numel, conn.connectome.synIdx(excSomaMask)), ...
       [height(axonMeta), 1]);
    axonMeta.excSomaSynFrac = ...
        axonMeta.excSomaSynCount ...
     ./ axonMeta.synCount;
    
    shaftPrefIds = find( ...
        axonMeta.synCount >= opts.minSynPre ...
      & axonMeta.fullPriSpineSynFrac < 0.5);
    unassignedIds = find( ...
        axonMeta.synCount >= opts.minSynPre & ~( ...
        axonMeta.fullPriSpineSynFrac >= 0.50 ...
      | axonMeta.fullPriSpineSynFrac <  0.25 ...
      | axonMeta.fullSecSpineSynFrac >  0.25 ...
      | axonMeta.excSomaSynFrac      >  0.10));
    tcCandIds = find( ...
        conn.axonMeta.fullPriSpineSynDens ...
      > opts.minSpineSynDensTc);
    
    % Sanity check
    assert(all(ismember(unassignedIds, shaftPrefIds)));

    % Excitatory axons
    axonClasses = struct;
    axonClasses(1).axonIds = find( ...
        axonMeta.synCount >= opts.minSynPre ...
      & axonMeta.fullPriSpineSynFrac >= 0.5);
    axonClasses(1).nullAxonIds = ...
        axonClasses(end).axonIds;
    axonClasses(1).title = sprintf( ...
        'excitatory axons with ≥ %d synapses (n = %d)', ...
        opts.minSynPre, numel(axonClasses(end).axonIds));

    % Inhibitory axons
    axonClasses(2).axonIds = setdiff( ...
        shaftPrefIds, unassignedIds);
    axonClasses(2).nullAxonIds = ...
        axonClasses(end).axonIds;
    axonClasses(2).title = sprintf( ...
        'inhibitory axons with ≥ %d synapses (n = %d)', ...
        opts.minSynPre, numel(axonClasses(end).axonIds));
    
    % Thalamocortical axons
    axonClasses(3).axonIds = intersect( ...
        axonClasses(1).axonIds, tcCandIds);
    axonClasses(3).nullAxonIds = axonClasses(1).axonIds;
    axonClasses(3).title = sprintf( ...
        'thalamocortical axons (n = %d)', ...
        numel(axonClasses(end).axonIds));
    
    % Corticocortical axons
    axonClasses(4).axonIds = setdiff( ...
        axonClasses(1).axonIds, tcCandIds);
    axonClasses(4).nullAxonIds = axonClasses(1).axonIds;
    axonClasses(4).title = sprintf( ...
        'corticocortical axons (n = %d)', ...
        numel(axonClasses(end).axonIds));
    
    % Update axon meta data
    conn.axonMeta.axonClass(:) = {'Other'};
    conn.axonMeta.axonClass(axonClasses(4).axonIds) = {'Corticocortical'};
    conn.axonMeta.axonClass(axonClasses(3).axonIds) = {'Thalamocortical'};
    conn.axonMeta.axonClass(axonClasses(2).axonIds) = {'Inhibitory'};
    
    axonClassOrder = { ...
        'Corticocortical', 'Thalamocortical', 'Inhibitory', 'Other'};
    conn.axonMeta.axonClass = categorical( ...
        conn.axonMeta.axonClass, axonClassOrder, 'Ordinal', true);
end
