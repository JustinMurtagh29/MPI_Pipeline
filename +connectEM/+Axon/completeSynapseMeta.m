function axonMeta = completeSynapseMeta( ...
        param, interSyn, boutonMeta, conn, syn)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    import connectEM.Axon.*;
    
    axonMeta = conn.axonMeta;
    synapses = syn.synapses;
    
    synapses.id = reshape( ...
        1:size(synapses, 1), [], 1);
    synapses.ontoSpine = ...
        synapses.type == 'PrimarySpine' ...
      | synapses.type == 'SecondarySpine';
  
    maxSegId = Seg.Global.getMaxSegId(param);
    axonLUT = Agglo.buildLUT(maxSegId, conn.axons);
    
    synapses.axonId = cellfun( ...
        @(segIds) setdiff(axonLUT(segIds), 0), ...
        synapses.presynId, 'UniformOutput', false);
    
    synapses(~cellfun(@isscalar, synapses.axonId), :) = [];
    synapses.axonId = cell2mat(synapses.axonId);
    
    if ~isempty(interSyn)
        axonMeta.pathLen = interSyn.axonPathLens / 1E3;
    else
        axonMeta.pathLen(:) = nan;
    end
    
    axonMeta.fullSynCount = accumarray( ...
        synapses.axonId, 1, size(axonMeta.id));
    axonMeta.fullSpineSynCount = accumarray( ...
        synapses.axonId, synapses.ontoSpine, size(axonMeta.id));
    
    priSynSpineMask = ...
        synapses.type == 'PrimarySpine';
    axonMeta.fullPriSpineSynCount = accumarray( ...
        synapses.axonId, priSynSpineMask, size(axonMeta.id));
    
    if ~isempty(interSyn)
        synIds = getSynapses(conn, syn);
        boutonIds = clusterSynapsesIntoBoutons(synIds, interSyn);
        fullPriSynSpineMask = syn.synapses.type == 'PrimarySpine';
        
        % Average number of primary spine innervations per axonal bouton
        axonMeta.fullPriSpinesPerBouton = cellfun( ...
            @(synIds, boutonIds) mean(accumarray( ...
                boutonIds, fullPriSynSpineMask(synIds))), ...
            synIds, boutonIds);

        % Fraction of axonal boutons with multiple primary spine synapses
        axonMeta.fullPriSpinesMultiHitFrac = cellfun( ...
            @(synIds, boutonIds) mean(accumarray( ...
                boutonIds, fullPriSynSpineMask(synIds)) > 1), ...
            synIds, boutonIds);
    else
        axonMeta.fullPriSpinesPerBouton(:) = nan;
        axonMeta.fullPriSpinesMultiHitFrac(:) = nan;
    end
    
    if ~isempty(boutonMeta)
        % Median volume of axonal boutons (in µm³)
        axonMeta.medianBoutonVol = cellfun( ...
            @median, boutonMeta.boutonVols) / (1E3) ^ 3;
    else
        axonMeta.medianBoutonVol(:) = nan;
    end
end
