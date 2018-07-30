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
 
    tcCandIds = find( ...
        conn.axonMeta.fullPriSpineSynDens ...
      > opts.minSpineSynDensTc);
  
    % For the source of these logistic regression parameters, see
    % connectEM.Axon.Script.detectThalamocorticals
    % 1f5d888c36ed0776851732be1cf36aafa6953b3b
    tcLogRegFeatNames = { ...
        'fullPriSpineSynDens', 'fullPriSpinesPerBouton', ...
        'fullPriSpinesMultiHitFrac', 'medianBoutonVol'};
    tcLogRegWeights = [0.356965; 1.664203; 0.725779; 0.829174];
    tcLogRegBias = -3.093222;
    
    tcProb = table2array(conn.axonMeta(:, tcLogRegFeatNames));
    tcProb = tcProb * tcLogRegWeights + tcLogRegBias;
    tcProb = 1 ./ (1 + exp(-tcProb));
    conn.axonMeta.tcProb = tcProb;

    % Excitatory axons
    axonClasses = struct;
    axonClasses(1).axonIds = find( ...
        conn.axonMeta.synCount >= opts.minSynPre ...
      & conn.axonMeta.fullPriSpineSynFrac >= 0.5);
    axonClasses(1).nullAxonIds = axonClasses(end).axonIds;
    axonClasses(1).title = sprintf( ...
        'excitatory axons with ≥ %d synapses (n = %d)', ...
        opts.minSynPre, numel(axonClasses(end).axonIds));

    % Inhibitory axons
    axonClasses(2).axonIds = find( ...
        conn.axonMeta.synCount >= opts.minSynPre ...
      & conn.axonMeta.fullPriSpineSynFrac < 0.2);
    axonClasses(2).nullAxonIds = axonClasses(end).axonIds;
    axonClasses(2).title = sprintf( ...
        'inhibitory axons with ≥ %d synapses (n = %d)', ...
        opts.minSynPre, numel(axonClasses(end).axonIds));
    
    % Thalamocortical axons
    axonClasses(3).axonIds = intersect( ...
        axonClasses(1).axonIds, tcCandIds);
    axonClasses(3).nullAxonIds = axonClasses(end).axonIds;
    axonClasses(3).title = sprintf( ...
        'thalamocortical axons (n = %d)', ...
        numel(axonClasses(end).axonIds));
    
    % Corticocortical axons
    axonClasses(4).axonIds = setdiff( ...
        axonClasses(1).axonIds, tcCandIds);
    axonClasses(4).nullAxonIds = axonClasses(end).axonIds;
    axonClasses(4).title = sprintf( ...
        'corticocortical axons (n = %d)', ...
        numel(axonClasses(end).axonIds));

    % Unclear axons
    axonClasses(5).axonIds = find( ...
        conn.axonMeta.synCount >= opts.minSynPre ...
      & conn.axonMeta.fullPriSpineSynFrac >= 0.2 ...
      & conn.axonMeta.fullPriSpineSynFrac  < 0.5);
    axonClasses(5).nullAxonIds = axonClasses(end).axonIds;
    axonClasses(5).title = sprintf( ...
        'unclear axons with ≥ %d synapses (n = %d)', ...
        opts.minSynPre, numel(axonClasses(end).axonIds));
    
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
