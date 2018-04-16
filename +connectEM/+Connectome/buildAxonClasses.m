function axonClasses = buildAxonClasses(conn, varargin)
    % axonClasses = buildAxonClasses(conn, varargin)
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opts = struct;
    opts.minSynPre = 1;
    opts = Util.modifyStruct(opts, varargin{:});
    
    conn.axonMeta.fullPriSpineSynFrac = ...
        conn.axonMeta.fullPriSpineSynCount ...
     ./ conn.axonMeta.fullSynCount;

    % Excitatory axons
    axonClasses = struct;
    axonClasses(1).axonIds = find( ...
        conn.axonMeta.synCount >= opts.minSynPre ...
      & conn.axonMeta.fullPriSpineSynFrac > 0.5);
    axonClasses(1).nullAxonIds = ...
        axonClasses(end).axonIds;
    axonClasses(1).title = sprintf( ...
        'excitatory axons with ≥ %d synapses (n = %d)', ...
        opts.minSynPre, numel(axonClasses(end).axonIds));
    
    axonClasses(1).specs = struct;
    axonClasses(1).specs.WholeCell.pThresh = 0.05;
    axonClasses(1).specs.ApicalDendrite.pThresh = 0.04;

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
    axonClasses(2).specs.Somata.pThresh = 0.11;
    axonClasses(2).specs.WholeCell.pThresh = 0.12;
    axonClasses(2).specs.ApicalDendrite.pThresh = 0.05;
    axonClasses(2).specs.SmoothDendrite.pThresh = 0.03;
    
    % Thalamocortical axons
    axonClasses(3).axonIds = find( ...
        conn.axonMeta.isThalamocortical);
    axonClasses(3).nullAxonIds = find( ...
        conn.axonMeta.synCount >= opts.minSynPre ...
      & conn.axonMeta.fullPriSpineSynFrac > 0.5);
    axonClasses(3).title = sprintf( ...
        'thalamocortical axons (n = %d)', ...
        numel(axonClasses(end).axonIds));
    
    axonClasses(3).specs = struct;
    axonClasses(3).specs.Somata.pThresh = 0.015;
    axonClasses(3).specs.WholeCell.pThresh = 0.195;
    axonClasses(3).specs.ApicalDendrite.pThresh = 0.025;
    
    % Corticocortical axons
    axonClasses(4).axonIds = setdiff( ...
        axonClasses(1).axonIds, ... % excitatory axons
        axonClasses(3).axonIds); % thalamocortical axons
    axonClasses(4).nullAxonIds = find( ...
        conn.axonMeta.synCount >= opts.minSynPre ...
      & conn.axonMeta.fullPriSpineSynFrac > 0.5);
    axonClasses(4).title = sprintf( ...
        'corticocortical axons (n = %d)', ...
        numel(axonClasses(end).axonIds));

    axonClasses(4).specs = struct;
    axonClasses(4).specs.WholeCell.pThresh = 0.04;
    axonClasses(4).specs.ApicalDendrite.pThresh = 0.04;
end
