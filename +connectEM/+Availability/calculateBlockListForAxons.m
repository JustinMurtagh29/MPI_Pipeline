function axonBlocks = calculateBlockListForAxons(conn, aggloSurfAreas)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    aggloIds = cellfun( ...
        @(asa) reshape(asa(:, 1:2), [], 1), ...
        aggloSurfAreas(:), 'UniformOutput', false);
    
    blockIds = repelem( ...
        reshape(1:numel(aggloIds), [], 1), ...
        cellfun(@(a) numel(a), aggloIds(:)));
    
    blockCoords = cell(1, 3);
   [blockCoords{:}] = ind2sub( ...
        size(aggloSurfAreas), blockIds);
    clear blockIds;
    
    t = table;
    t.aggloId = cell2mat(aggloIds);
    t.block = cell2mat(blockCoords);
    
    t(t.aggloId >= 0, :) = [];
    t.aggloId = -t.aggloId;
    t = unique(t, 'rows');
    
    axonBlocks = accumarray( ...
        t.aggloId, (1:size(t, 1))', size(conn.axons), ...
        @(rows) {t.block(rows, :)}, {zeros(0, 3)});
end