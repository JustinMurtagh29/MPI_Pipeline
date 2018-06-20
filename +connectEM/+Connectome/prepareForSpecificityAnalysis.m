function [conn, axonClasses] = ...
        prepareForSpecificityAnalysis(conn, axonClasses, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opt = struct;
    opt.minSynPre = 1;
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
    
    %% Fix order of categories
    cats = unique(conn.denMeta.targetClass);
    conn.denMeta.targetClass = categorical( ...
        conn.denMeta.targetClass, cats);
    conn.denMeta.targetClass = reordercats( ...
        conn.denMeta.targetClass, targetClassOrder);
    
    %% Remove excitatory synapses onto excitatory somata
    % As discussed with MH in a meeting on 19.06.2018.
    conn = removeExcSynapsesOntoExcSomata(conn);
    
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
    
    % Decrease synapse count
    synT = table;
    synT.id = cell2mat(excSomaConn.synIdx);
    synT.axonId(:) = repelem( ...
        excSomaConn.edges(:, 1), ...
        cellfun(@numel, excSomaConn.synIdx));
    synT = unique(synT, 'rows');
    
   [excAxonRows, ~, excAxonSynCount] = unique(synT.axonId);
   [~, excAxonRows] = ismember(excAxonRows, conn.axonMeta.id);
    excAxonSynCount = accumarray(excAxonSynCount, 1);
    
    conn.axonMeta.synCount(excAxonRows) = ...
        conn.axonMeta.synCount(excAxonRows) - excAxonSynCount;
    
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
    
    % Same for connectome meta data
    conn = rmfield(conn, 'connectomeMeta');
end
