function pairConfigs = buildPairConfigs(synT, plotConfig)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    synT.id = reshape(1:size(synT, 1), [], 1);
    synT = synT(plotConfig.synIds, :);
    
   [~, uniSynT, synT.pairId] = unique( ...
        synT(:, {'preAggloId', 'postAggloId'}), 'rows');
    couplings = accumarray(synT.pairId, 1);
    uniSynT = synT(uniSynT, :);
    
    % Same-axon same-dendrite
    sameSameIds = synT;
    sameSameIds = sortrows(sameSameIds, 'pairId');
    sameSameIds.coupling = couplings(sameSameIds.pairId);
    sameSameIds(sameSameIds.coupling ~= 2, :) = [];
    sameSameIds = reshape(sameSameIds.id, 2, [])';
    
    % Same-axon different-dendrite
    sameDiffIds = cell2mat(accumarray( ...
        uniSynT.preAggloId, uniSynT.id, [], @(ids) { ...
        reshape(ids(1:(2 * floor(numel(ids) / 2))), [], 2)}));
    
    % Different-axon same-dendrite
    diffSameIds = cell2mat(accumarray( ...
        uniSynT.postAggloId, uniSynT.id, [], @(ids) { ...
        reshape(ids(1:(2 * floor(numel(ids) / 2))), [], 2)}));
    
    % Different-axon different-dendrite
    rng(0);
    diffDiffIds = randperm(2 * floor(height(synT) / 2));
    diffDiffIds = reshape(diffDiffIds, [], 2);
    
    diffDiffIds = diffDiffIds( ...
        diff(synT.preAggloId(diffDiffIds), 1, 2) ...
      & diff(synT.postAggloId(diffDiffIds), 1, 2), :);
    diffDiffIds = synT.id(diffDiffIds);
    
    % Random pairs
    rng(0);
    randPairIds = randperm(2 * floor(height(synT) / 2));
    randPairIds = synT.id(reshape(randPairIds, [], 2));
    
    % Build output
    pairConfigs = struct;
    pairConfigs(1).synIdPairs = sameSameIds;
    pairConfigs(1).title = 'Same-axon same-dendrite';
    
    pairConfigs(2).synIdPairs = sameDiffIds;
    pairConfigs(2).title = 'Same-axon different-dendrite';
    
    pairConfigs(3).synIdPairs = diffSameIds;
    pairConfigs(3).title = 'Different-axon same-dendrite';
    
    pairConfigs(4).synIdPairs = diffDiffIds;
    pairConfigs(4).title = 'Different-axon different-dendrite';
    
    pairConfigs(5).synIdPairs = randPairIds;
    pairConfigs(5).title = 'Random pairs';
    
    % Sort synapse by size
    synIdPairs = arrayfun( ...
        @(p) sortByArea(synT, p.synIdPairs), ...
        pairConfigs, 'UniformOutput', false);
   [pairConfigs.synIdPairs] = deal(synIdPairs{:});
    
    titles = arrayfun(@(c) sprintf( ...
        '%s (n = %d)', c.title, size(c.synIdPairs, 1)), ...
        pairConfigs, 'UniformOutput', false);
   [pairConfigs.title] = deal(titles{:});
end

function idPairs = sortByArea(synT, idPairs)
    idPairs = transpose(idPairs);
   [~, sortIds] = ismember(idPairs, synT.id);
   [~, sortIds] = sort(synT.area(sortIds), 1, 'descend');
    sortIds = sortIds + 2 * ((1:size(sortIds, 2)) - 1);
    idPairs = transpose(idPairs(sortIds));
end
