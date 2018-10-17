function [conn, axonClasses] = ...
        prepareForSpecificityAnalysis(conn, axonClasses, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opt = struct;
    opt.syn = [];
    opt.minSynPre = 1;
    opt.perisomaClass = false;
    opt.removeSecondarySpines = false;
    opt = Util.modifyStruct(opt, varargin{:});
    
    %% Configuration
    targetClassOrder = { ...
        'Somata', ...
        'ProximalDendrite', ...
        'SmoothDendrite', ...
        'ApicalDendrite', ...
        'AxonInitialSegment', ...
        'OtherDendrite'};
    
    %% Treat soma-based reconstructions
    inMask = conn.denMeta.isInterneuron;
    wcMask = conn.denMeta.targetClass == 'WholeCell';

    % Inhibitory whole cell → smooth dendrite
    conn.denMeta.targetClass(wcMask &  inMask) = 'SmoothDendrite';
    
    % Excitatory whole cell → proximal dendrite
    conn.denMeta.targetClass(wcMask & ~inMask) = 'ProximalDendrite';
    
    %% Remove excitatory synapses onto excitatory somata
    % As discussed with MH in a meeting on 19.06.2018.
    conn = removeExcSynapsesOntoExcSomata(conn);
    
    %% Lump somata and proximal dendrites into a "perisomatic" class
    % This option was added as a sanity for the specificity analysis. The
    % question is aims to address is the following:
    %   How often does it happen that an inhibitory axon has hints of
    % specificity for both somata and proximal dendrites, but not enough
    % to be detected as specific for any of those classes.
    %   Here, we try to answer this question by lumping both somata and
    % proximal dendrites into a single class, called 'Perisoma'.
    if opt.perisomaClass
        perisomaClasses = {'Somata', 'ProximalDendrite'};
        
        curMask = ismember(conn.denMeta.targetClass, perisomaClasses);
        conn.denMeta.targetClass(curMask) = {'Perisoma'};
        
        targetClassOrder(ismember(targetClassOrder, perisomaClasses)) = [];
        targetClassOrder = horzcat({'Perisoma'}, targetClassOrder);
    end
    
    %% Remove secondary spine innervations
    % As discussed with MH in a meeting on 05.07.2018.
    if opt.removeSecondarySpines
        assert(~isempty(opt.syn));
        conn = removeSecondarySpineSynapses(conn, opt.syn);
    end
    
    %% Update meta data (after potential synapse removal)
    conn = updateMetaData(conn);
    
    %% Fix order of categories
    cats = unique(conn.denMeta.targetClass);
    conn.denMeta.targetClass = categorical( ...
        conn.denMeta.targetClass, cats);
    conn.denMeta.targetClass = reordercats( ...
        conn.denMeta.targetClass, targetClassOrder);
    
    %% Update axon classes, if desired
    if ~exist('axonClasses', 'var'); return; end
    
    validAxonIds = conn.axonMeta.synCount >= opt.minSynPre;
    validAxonIds = conn.axonMeta.id(validAxonIds);
    
    for curIdx = 1:numel(axonClasses)
        axonClasses(curIdx).axonIds = intersect( ...
            axonClasses(curIdx).axonIds, validAxonIds);
        axonClasses(curIdx).nullAxonIds = intersect( ...
            axonClasses(curIdx).nullAxonIds, validAxonIds);
        axonClasses(curIdx).title = regexprep( ...
            axonClasses(curIdx).title, 'n = \d+', ...
            sprintf('n = %d', numel(axonClasses(curIdx).axonIds)));
    end
end

function conn = removeSecondarySpineSynapses(conn, syn)
    synMask = syn.synapses.type ~= 'SecondarySpine';
    
    conn.connectome.synIdx = cellfun( ...
        @(synIds) synIds(synMask(synIds)), ...
        conn.connectome.synIdx, 'UniformOutput', false);
    conn.connectome(cellfun( ...
        @isempty, conn.connectome.synIdx), :) = [];
end

function conn = removeExcSynapsesOntoExcSomata(conn)
    excAxonIds = ismember( ...
        conn.axonMeta.axonClass, ...
        {'Corticocortical', 'Thalamocortical'});
    excAxonIds = conn.axonMeta.id(excAxonIds);
    
    excSomaIds = ...
        ~conn.denMeta.isInterneuron ...
       & conn.denMeta.targetClass == 'Somata';
    excSomaIds = conn.denMeta.id(excSomaIds);
    
    excSomaConn = conn.connectome;
    excSomaConn.row(:) = 1:height(excSomaConn);
    excSomaConn = excSomaConn(...
        ismember(excSomaConn.edges(:, 1), excAxonIds) ...
      & ismember(excSomaConn.edges(:, 2), excSomaIds), :);
    
    % Remove synapses from connectome
    conn.connectome(excSomaConn.row, :) = [];
end

function conn = updateMetaData(conn)
    axonIds = repelem( ...
        conn.connectome.edges(:, 1), ...
        cellfun(@numel, conn.connectome.synIdx));
   [~, axonIds] = ismember(axonIds, conn.axonMeta.id);
    conn.axonMeta.synCount = accumarray( ...
        axonIds, 1, [numel(conn.axons), 1]);
    
    % HACK(amotta): I'm too lazy to update the `spineSynCount`,
    % `fullSynCount`, `fullSpineSynCount`, and `fullPriSpineSynCount`
    % variables. To make sure that nobody accidentally uses wrong numbers,
    % let's just remove the variables completely.
    conn.axonMeta(:, { ...
        'spineSynCount', ...
        'fullSynCount', ...
        'fullSpineSynCount', ...
        'fullPriSpineSynCount', ...
        'fullPriSpineSynFrac', ...
        'fullPriSpineSynDens'}) = [];
    
    conn.denMeta(:, { ...
        'synCount', ...
        'spineSynCount'}) = [];
    
    % Same for connectome meta data
    conn = rmfield(conn, 'connectomeMeta');
end
