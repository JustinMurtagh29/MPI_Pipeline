function synTypes = classifyType( ...
        param, synapses, synScores, shAgglos, somaAgglos)
    % synTypes = classifyType( ...
    %     param, synapses, synProbs, shAgglos, somaAgglos)
    %
    %   Classifies synapses into
    %   * shaft synapses
    %   * soma synapses
    %   * primary spine innervations
    %   * secondary spine innervations
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
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

    % Sanity check
    % Synapse threshold was -1.67
    assert(all(cellfun(@min, synapses.synScores) > -1.67));
    synapses.maxSynScore = cellfun(@max, synapses.synScores);

    synapses.shId = cellfun( ...
        @(segIds) max(shLUT(segIds)), ...
        synapses.postsynId);
    synapses.somaId = cellfun( ...
        @(segIds) max(somaLUT(segIds)), ...
        synapses.postsynId);

    %% Distringuish primary and secondary spine synapses
    shSynIds = accumarray( ...
        1 + synapses.shId, synapses.id, ...
       [1 + numel(shAgglos), 1], @(ids) {ids}, {zeros(0, 1)});
    shSynIds = shSynIds(2:end);

    % Sort synapse by SynEM scores
    % Most probable synapse becomes primary
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
    
    %% Build output
    synTypes = cell(size(synapses.id));
    synTypes(:) = {'Shaft'};
    synTypes(shT.priSynId) = {'PrimarySpine'};
    synTypes(cell2mat(shT.secSynIds)) = {'SecondarySpine'};
    synTypes(synapses.somaId ~= 0) = {'Soma'};
    synTypes = categorical(synTypes);
end