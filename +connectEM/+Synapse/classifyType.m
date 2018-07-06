function synTypes = classifyType( ...
        param, synapses, synScores, shAgglos, somaAgglos, axonAgglos)
    % synTypes = classifyType( ...
    %     param, synapses, synProbs, shAgglos, somaAgglos, axonAgglos)
    %
    %   Classifies synapses into
    %   * shaft synapses
    %   * soma synapses
    %   * primary spine innervations
    %   * secondary spine innervations
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    if ~exist('axonAgglos', 'var')
        axonAgglos = cell(0, 1);
    end
    
    %% Augment synapse table
    maxSegId = Seg.Global.getMaxSegId(param);
    
    % Building look-up tables
    % Soma class dominates over spine heads
    somaLUT = Agglo.buildLUT(maxSegId, somaAgglos);
    shLUT = Agglo.buildLUT(maxSegId, shAgglos);
    shLUT(somaLUT ~= 0) = 0;

    synapses.id = reshape( ...
        1:size(synapses, 1), [], 1);
    synapses.synScores = cellfun( ...
        @(ids) max(synScores(ids, :), [], 2), ...
        synapses.edgeIdx, 'UniformOutput', false);
    synapses.maxSynScore = cellfun(@max, synapses.synScores);

    synapses.shId = cellfun( ...
        @(segIds) max(shLUT(segIds)), ...
        synapses.postsynId);
    synapses.somaId = cellfun( ...
        @(segIds) max(somaLUT(segIds)), ...
        synapses.postsynId);

    %% Distringuish primary and secondary spine synapses
    % For each spine head which is postsynaptic to at least one synapse, we
    % consider the synapse with the highest SynEM score to be the primary
    % spine innervation. All other synapses onto that spine head are
    % secondary innervations.
    
    shSynIds = accumarray( ...
        1 + synapses.shId, synapses.id, ...
       [1 + numel(shAgglos), 1], @(ids) {ids}, {zeros(0, 1)});
    shSynIds = shSynIds(2:end);

    % Sort synapse by SynEM scores
   [~, sortedIds] = cellfun( ...
        @(synIds) sort(synapses.maxSynScore(synIds)), ...
        shSynIds, 'UniformOutput', false);
    shSynIds = cellfun( ...
        @(synIds, sortIds) synIds(sortIds), ...
        shSynIds, sortedIds, 'UniformOutput', false);
    clear sortedIds;

    shT = table;
    shT.id = reshape( ...
        1:numel(shAgglos), [], 1);

    shT(cellfun(@isempty, shSynIds), :) = [];
    shSynIds(cellfun(@isempty, shSynIds)) = [];

    shT.priSynId = cellfun( ...
        @(synIds) synIds(end), shSynIds);
    shT.secSynIds = cellfun( ...
        @(synIds) reshape(synIds(1:(end - 1)), [], 1), ...
        shSynIds, 'UniformOutput', false);
    
    priSpineSynIds = shT.priSynId;
    secSpineSynIds = cell2mat(shT.secSynIds);
    
    %% If possible, handle split axon-spine interfaces
    % If axon A has a primary spine synapse onto spine head SH, then all
    % synapses from A onto SH are marked as primary spine synapses.
    
    if ~isempty(axonAgglos)
        axonLUT = Agglo.buildLUT(maxSegId, axonAgglos);
        
        synapses.axonId = cellfun( ...
            @(segIds) max(axonLUT(segIds)), ...
            synapses.presynId);
        synapses.axonId(~synapses.axonId) = ...
            numel(axonAgglos) + (1:sum(~synapses.axonId));
        
        priSpineSynIds = synapses(ismember( ...
            synapses.id, priSpineSynIds), :);
        priSpineSynIds = synapses.id(ismember( ...
            synapses(:,  {'shId', 'axonId'}), ...
            priSpineSynIds(:, {'shId', 'axonId'}), 'rows'));
        secSpineSynIds = setdiff(secSpineSynIds, priSpineSynIds);
    end
    
    % Sanity check
    assert(isempty(intersect( ...
        priSpineSynIds, secSpineSynIds)));
    
    %% Build output
    synTypes = cell(size(synapses, 1), 1);
    synTypes(:) = {'Shaft'};
    synTypes(priSpineSynIds) = {'PrimarySpine'};
    synTypes(secSpineSynIds) = {'SecondarySpine'};
    synTypes(synapses.somaId ~= 0) = {'Soma'};
    synTypes = categorical(synTypes);
end
