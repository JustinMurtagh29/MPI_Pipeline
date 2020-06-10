function [connLin, synLin, connLinFile] = ...
        loadConnectome(param, connFile, connLinFile, classDef)
    % loadConnectome(param)
    %   This is a wrapper around the standard connectEM.Connectome.load
    %   function specifically for the analysis of synapse size consistency.
    %   For the consistency analysis, we use axon reconstructions that have
    %   been linearized by splitting up all branch points. This reduces the
    %   risk of wrong same-axon same-dendrite synapses caused by mergers.
    %
    %   But the linearization renders axon class inference more difficult.
    %   Hence, we want to reuse the axon classes from before splitting.
    %   This function loads the connectome based on fully linearized axons,
    %   but re-uses the previously determined axon classes.
    %
    % Written by
    %   Alessandro motta <alessandro.motta@brain.mpg.de>
    if ~exist('classDef', 'var') || isempty(classDef)
        classDef = 'axonClasses';
    end

    %% Loading data
   [conn, ~, axonClasses] = connectEM.Connectome.load(param, connFile);
   [connLin, synLin] = connectEM.Connectome.load(param, connLinFile);
    
    % Add parent IDs for inheritance of axon classes
    curMeta = load(conn.info.param.axonFile);
    conn.axonMeta.unsplitParentId = ...
        curMeta.parentIds(conn.axonMeta.parentId);
    
    curMeta = load(connLin.info.param.axonFile);
    connLin.axonMeta.unsplitParentId = ...
        curMeta.parentIds(connLin.axonMeta.parentId);
    
    %% Transfer axon classes
    connLin = ...
        connectEM.Consistency.loadConnectomeCore( ...
            classDef, conn, axonClasses, connLin);
end
