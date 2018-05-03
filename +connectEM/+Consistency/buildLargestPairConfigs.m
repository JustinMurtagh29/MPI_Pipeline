function pairConfigs = ...
        buildLargestPairConfigs(synT, plotConfig, coupling, nthLargest)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    if ~exist('nthLargest', 'var') || isempty(nthLargest)
        nthLargest = 1;
    end
    
    synT.id = reshape(1:height(synT), [], 1);
    synT = synT(plotConfig.synIds, :);
    
    % Control
    rng(0);
    ctrlIds = coupling * floor(height(synT) / coupling);
    ctrlIds = reshape(randperm(ctrlIds), coupling, []);
    
   [~, sortIds] = sort(synT.area(ctrlIds), 'descend');
    sortIds = sortIds + (1:coupling:numel(sortIds)) - 1;
    sortIds = sortIds(nthLargest:(nthLargest + 1), :);
    ctrlIds = transpose(synT.id(ctrlIds(sortIds)));
    
    % Observed
   [~, ~, synT.pairId] = unique(synT(:, ...
        {'preAggloId', 'postAggloId'}), 'rows');
    couplings = accumarray(synT.pairId, 1);
    synT.coupling = couplings(synT.pairId);
    
    synT(synT.coupling ~= coupling, :) = [];
    synT = sortrows(synT, {'pairId', 'area'}, 'descend');
    
    synIds = reshape(synT.id, coupling, []);
    synIds = synIds(nthLargest:(nthLargest + 1), :);
    synIds = transpose(synIds);
    
    % Build output
    pairConfigs = struct;
    pairConfigs(1).synIdPairs = synIds;
    pairConfigs(1).title = sprintf( ...
        'Observed %dth largest pairs (n = %d)', ...
        nthLargest, size(synIds, 1));
    pairConfigs(2).synIdPairs = ctrlIds;
    pairConfigs(2).title = sprintf( ...
        'Expected %dth largest pair (n = %d)', ...
        nthLargest, size(ctrlIds, 1));
end
