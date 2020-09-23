function [connLin, synLin, connLinFile] = ...
        loadConnectomePaper(param, classDef)
    % loadConnectomePaper(param)
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
    
    rootDir = param.saveFolder;
    
    % Connectome based on partially split axons, which were used to develop
    % and validate the criteria for the identification of TC axons.
    connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
    
    % Connectome based on fully linearized axons. These split axons reduce
    % the risk of wrong "same-axon same-dendrite" synapses introduced by
    % merge errors.
    connLinFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

    %% Loading data
    param = load(fullfile(rootDir, 'allParameter.mat'));
    param = param.p;

   [conn, ~, axonClasses] = connectEM.Connectome.load(param, connFile);
   [connLin, synLin] = connectEM.Connectome.load(param, connLinFile);

    % HACK(amotta): The connectome used `axons_19_a_partiallySplit_v2.mat`,
    % but the parent IDs were only added later in version 4...
    curMeta = conn.info.param.axonFile;
    curMeta = load(strrep(curMeta, '_v2.mat', '_v4.mat'));
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
