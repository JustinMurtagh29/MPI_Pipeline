function axonClasses = buildAxonClasses(conn, varargin)
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
 
    tcCandIds = find( ...
        conn.axonMeta.fullPriSpineSynDens ...
      > opts.minSpineSynDensTc);

    % Excitatory axons
    axonClasses = struct;
    axonClasses(1).axonIds = find( ...
        conn.axonMeta.synCount >= opts.minSynPre ...
      & conn.axonMeta.fullPriSpineSynFrac >= 0.5);
    axonClasses(1).nullAxonIds = ...
        axonClasses(end).axonIds;
    axonClasses(1).title = sprintf( ...
        'excitatory axons with ≥ %d synapses (n = %d)', ...
        opts.minSynPre, numel(axonClasses(end).axonIds));
    
    axonClasses(1).specs = struct;
    axonClasses(1).specs.WholeCell.pThresh = 0.05;
    axonClasses(1).specs.ApicalDendrite.pThresh = 0.035;

    % Inhibitory axons
    axonClasses(2).axonIds = find( ...
        conn.axonMeta.synCount >= opts.minSynPre ...
      & conn.axonMeta.fullPriSpineSynFrac < 0.5);
    axonClasses(2).nullAxonIds = ...
        axonClasses(end).axonIds;
    axonClasses(2).title = sprintf( ...
        'inhibitory axons with ≥ %d synapses (n = %d)', ...
        opts.minSynPre, numel(axonClasses(end).axonIds));
    
    axonClasses(2).specs = struct;
    axonClasses(2).specs.Somata.pThresh = 0.115;
    axonClasses(2).specs.WholeCell.pThresh = 0.12;
    axonClasses(2).specs.ApicalDendrite.pThresh = 0.06;
    axonClasses(2).specs.SmoothDendrite.pThresh = 0.055;
    
    % Thalamocortical axons
    axonClasses(3).axonIds = intersect( ...
        axonClasses(1).axonIds, tcCandIds);
    axonClasses(3).nullAxonIds = axonClasses(1).axonIds;
    axonClasses(3).title = sprintf( ...
        'thalamocortical axons (n = %d)', ...
        numel(axonClasses(end).axonIds));
    
    axonClasses(3).specs = struct;
    axonClasses(3).specs.WholeCell.pThresh = 0.215;
    axonClasses(3).specs.ApicalDendrite.pThresh = 0.025;
    
    % Corticocortical axons
    axonClasses(4).axonIds = setdiff( ...
        axonClasses(1).axonIds, tcCandIds);
    axonClasses(4).nullAxonIds = axonClasses(1).axonIds;
    axonClasses(4).title = sprintf( ...
        'corticocortical axons (n = %d)', ...
        numel(axonClasses(end).axonIds));

    axonClasses(4).specs = struct;
    axonClasses(4).specs.WholeCell.pThresh = 0.04;
    axonClasses(4).specs.ApicalDendrite.pThresh = 0.04;
end
