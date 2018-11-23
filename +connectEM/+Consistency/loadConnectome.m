function [connLin, synLin, connLinFile] = loadConnectome(param, classDef)
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

    % HACK(amotta): The connectome used `axons_19_a_partiallySplit_v2.mat`,
    % but the parent IDs were only added later in version 4...
    curMeta = load(strrep(conn.info.param.axonFile, '_v2.mat', '_v4.mat'));
    conn.axonMeta.unsplitParentId = ...
        curMeta.parentIds(conn.axonMeta.parentId);

    %% Define axon classes
    switch classDef
        case 'axonClasses'
            conn.axonMeta.axonClass = renamecats( ...
                conn.axonMeta.axonClass, 'Inhibitory', 'Ignore');
            
        case 'specificityClasses'
           [conn, axonClasses] = ...
                connectEM.Connectome.prepareForSpecificityAnalysis( ...
                    conn, axonClasses, 'minSynPre', 10);
            curClasses = ...
                connectEM.Connectome.buildAxonSpecificityClasses( ...
                    conn, axonClasses(1:2));
            curClasses = curClasses(1).specs;

            curNames = fieldnames(curClasses);
            curAxonIds = cellfun( ...
                @(n) curClasses.(n).axonIds, ...
                curNames, 'UniformOutput', false);

            % Make specificity classes mutually non-overlapping
            curDupAxonIds = sort(cell2mat(curAxonIds));
            curDupAxonIds = curDupAxonIds(diff(curDupAxonIds) == 0);

            curAxonIds = cellfun( ...
                @(ids) setdiff(ids, curDupAxonIds), ...
                curAxonIds, 'UniformOutput', false);
            curExcAxonIds = axonClasses(1).axonIds;

            % NOTE(amotta): All non-excitatory axons (i.e., inhibitory and
            % unidentified axons) are marked "ignore". Excitatory axons
            % with no target specificity will be in the "other" category.
            conn.axonMeta.axonClass = repelem( ...
                {'Ignore'}, numel(conn.axons), 1);
            conn.axonMeta.axonClass(curExcAxonIds) = {'Other'};
            conn.axonMeta.axonClass(cell2mat(curAxonIds)) = ...
                repelem(curNames, cellfun(@numel, curAxonIds), 1);
            conn.axonMeta.axonClass = categorical(conn.axonMeta.axonClass);
            
        otherwise
            error('Unknown class definition "%s"', classDef);
    end
    
    %% Translate classes from partially split to fully linearized axons
    % NOTE(amotta): The partially split and the linearized version of the
    % axons were both derived from axons 19a. To translate the axon classes
    % from the partially split to the fully linearized axons, we thus have
    % to do a detour via the unsplit version.
    %   Because there's no bijective mapping between partially split and
    % fully linearized axons, we only look at the axons 19a which were
    % uniquely classified in the partially split version.
    clear cur*;
    
    curAllCats = conn.axonMeta.axonClass;
    curOtherCat = cat(1, curAllCats([]), {'Other'});
    curCats = setdiff(curAllCats, {'Other'; 'Ignore'});
    
    curParentIds = find(accumarray( ...
        conn.axonMeta.unsplitParentId, ...
        double(conn.axonMeta.axonClass), [], ...
        @(t) isscalar(setdiff(t, double(curOtherCat)))));
    curParentIds = arrayfun( ...
        @(c) intersect(conn.axonMeta.unsplitParentId( ...
            conn.axonMeta.axonClass == c), curParentIds), ...
        curCats, 'UniformOutput', false);
    
    % Make sure that parent lists are mutually non-overlapping
    Agglo.check(curParentIds);
    
   [connLin, synLin] = ...
       connectEM.Connectome.load(param, connLinFile);
   
    curMeta = load(connLin.info.param.axonFile);
    connLin.axonMeta.unsplitParentId = ...
        curMeta.parentIds(connLin.axonMeta.parentId);

    % Assign axon classes to linearized axons
    connLin.axonMeta.axonClass = repelem( ...
        curOtherCat, numel(connLin.axons), 1);
    for curIdx = 1:numel(curCats)
        connLin.axonMeta.axonClass(ismember( ...
            connLin.axonMeta.unsplitParentId, ...
            curParentIds{curIdx})) = curCats(curIdx);
    end
end
