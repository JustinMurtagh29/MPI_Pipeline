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
    
    % Random pairs
    rng(0);
    randIds = randperm(2 * floor(size(synT, 1) / 2));
    randIds = reshape(synT.id(randIds), [], 2);
    
    % Build output
    pairConfigs = struct;
    pairConfigs(1).synIdPairs = sameSameIds;
    pairConfigs(1).title = 'Same-axon same-dendrite';
    
    pairConfigs(2).synIdPairs = sameDiffIds;
    pairConfigs(2).title = 'Same-axon different-dendrite';
    
    pairConfigs(3).synIdPairs = diffSameIds;
    pairConfigs(3).title = 'Different-axon same-dendrite';
    
    pairConfigs(4).synIdPairs = randIds;
    pairConfigs(4).title = 'Random pairs';
    
    titles = arrayfun(@(c) sprintf( ...
        '%s (n = %d)', c.title, size(c.synIdPairs, 1)), ...
        pairConfigs, 'UniformOutput', false);
   [pairConfigs.title] = deal(titles{:});
end
