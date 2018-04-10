function axonMeta = completeSynapseMeta(param, conn, syn)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
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

    axonMeta = conn.axonMeta;
    axonMeta.fullSynCount = accumarray( ...
        synapses.axonId, 1, size(axonMeta.id));
    axonMeta.fullSpineSynCount = accumarray( ...
        synapses.axonId, synapses.ontoSpine, size(axonMeta.id));
    
    priSpineSynMask = ...
        synapses.type == 'PrimarySpine';
    axonMeta.fullPriSpineSynCount = accumarray( ...
        synapses.axonId, priSpineSynMask, size(axonMeta.id));
end