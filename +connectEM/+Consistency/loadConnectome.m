function [connLin, synLin, connLinFile] = loadConnectome(param)
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

    conn = connectEM.Connectome.load(param, connFile);

    % HACK(amotta): The connectome used `axons_19_a_partiallySplit_v2.mat`,
    % but the parent IDs were only added later in version 4...
    curMeta = load(strrep(conn.info.param.axonFile, '_v2.mat', '_v4.mat'));
    conn.axonMeta.unsplitParentId = ...
        curMeta.parentIds(conn.axonMeta.parentId);

   [connLin, synLin] = connectEM.Connectome.load(param, connLinFile);

    curMeta = load(connLin.info.param.axonFile);
    connLin.axonMeta.unsplitParentId = ...
        curMeta.parentIds(connLin.axonMeta.parentId);

    %% Translate classes from partially split to fully linearized axons
    % NOTE(amotta): The partially split and the linearized version of the
    % axons were both derived from axons 19a. To translate the axon classes
    % from the partially split to the fully linearized axons, we thus have
    % to do a detour via the unsplit version.
    %   Because there's no bijective mapping between partially split and
    % fully linearized axons, we only look at the axons 19a which were
    % uniquely classified in the partially split version.
    clear cur*;

    % Find parents of axons classified as TC / CC
    curTcParentIds = conn.axonMeta.unsplitParentId( ...
        conn.axonMeta.axonClass == 'Thalamocortical');
    curCcParentIds = conn.axonMeta.unsplitParentId( ...
        conn.axonMeta.axonClass == 'Corticocortical');

    % Check if there exist siblings of TC axons that were classified
    % differently (i.e., as corticocortical or inhibitory; other is okay).
    % If yes, we do not trust or use these TC axons.
    curTcSiblingsOkay = accumarray( ...
        conn.axonMeta.unsplitParentId, ...
        double(conn.axonMeta.axonClass), [], ...
        @(types) isequal(setdiff(types, 4), 2)); % 2 = TC, 4 = other
    curTcSiblingsOkay = curTcSiblingsOkay(curTcParentIds);
    curTcParentIds = curTcParentIds(curTcSiblingsOkay);

    curCcSiblingsOkay = accumarray( ...
        conn.axonMeta.unsplitParentId, ...
        double(conn.axonMeta.axonClass), [], ...
        @(types) isequal(setdiff(types, 4), 1)); % 1 = CC, 4 = other
    curCcSiblingsOkay = curCcSiblingsOkay(curCcParentIds);
    curCcParentIds = curCcParentIds(curCcSiblingsOkay);

    % Sanity check
    assert(isempty(intersect(curTcParentIds, curCcParentIds)));

    % Assign axon classes to linearized axons
    connLin.axonMeta.axonClass(:) = 'Other';
    connLin.axonMeta.axonClass(ismember( ...
        connLin.axonMeta.unsplitParentId, ...
        curTcParentIds)) = 'Thalamocortical';
    connLin.axonMeta.axonClass(ismember( ...
        connLin.axonMeta.unsplitParentId, ...
        curCcParentIds)) = 'Corticocortical';
    
    % Remove temporary field
    connLin.axonMeta(:, {'unsplitParentId'}) = [];
end
