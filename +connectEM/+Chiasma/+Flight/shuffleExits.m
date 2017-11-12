function [exits, rows] = shuffleExits(exits, batchSize)
    % [exits, rows] = shuffleExits(exits, batchSize)
    %   Shuffles the order of the exits to be queried. 
    %
    %   An optional `batchSize` can be specified. If this parameter is set,
    %   queries are shuffled in such a way that `batchSize` chiasmata are
    %   fully queried, before the next batch of chiasmata are processed.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    rng(0);
    exitsIn = exits;
    exits.row = reshape(1:size(exits, 1), [], 1);
    
    %% shuffle
    randIds = randperm(size(exits, 1));
    exits = exits(randIds, :);
    
    %% make batches
   [~, ~, uniChiasma] = unique( ...
       exits(:, {'aggloId', 'chiasmaId'}), 'rows', 'stable');
   
    if ~exist('batchSize', 'var') || isempty(batchSize)
        % default to no batching
        batchSize = max(uniChiasma);
    end
    
    for curIdx = batchSize:batchSize:max(uniChiasma)
        % select chiasmata below fold
        curMask = uniChiasma < curIdx;
        
        % bubble to the top
        exits = cat(1, exits(curMask, :), exits(~curMask, :));
        uniChiasma = cat(1, uniChiasma(curMask), uniChiasma(~curMask));
    end
    
    rows = exits.row;
    exits.row = [];
    
    %% sanity check
    assert(isequal(unique(exits, 'rows'), unique(exitsIn, 'rows')));
end