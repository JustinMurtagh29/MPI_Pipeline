function splitAgglos = splitAgglo( ...
        param, chiasmaParam, agglo, chiasmata, overlaps)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % sanity check for input arguments
    assert(numel(overlaps) == numel(chiasmata.ccCenterIdx));
    
    % build `queries` table for splitting
    queries = buildQueryTable(chiasmata, overlaps);
    queries(~queries.split, :) = [];
    queries.split = [];
    
    splitAgglos = ...
        connectEM.Chiasma.Ortho.splitAgglo( ...
            param, chiasmaParam, agglo, queries);
end

function queries = buildQueryTable(chiasmata, overlaps)
    % Builds `queries` table for the splitting routine
    queries = table;
    queries.centerNodeId = chiasmata.ccCenterIdx;
    queries.exitNodeIds = cellfun(@(ids) ...
        ids(:), chiasmata.queryIdx, 'UniformOutput', false);
    
   [queries.split, queries.exits] = cellfun( ...
       @buildExitTable, overlaps, 'UniformOutput', false);
    queries.split = cell2mat(queries.split);
end

function [split, exits] = buildExitTable(overlaps)
    nrExits = size(overlaps, 1);
    
    exits = table;
    exits.id(:) = (1:nrExits)';
    exits.groupId(:) = 0;
    
    % split only if all exits were visited at least once
    split = all(ismember(exits.id, overlaps));
    if ~split; return; end
    
    % build groups
    adjMat = sort(overlaps, 2);
    adjMat = adjMat(all(adjMat, 2), :);
    adjMat = sparse(adjMat(:, 2), adjMat(:, 1), 1, nrExits, nrExits);
    
   [~, groupIds] = graphconncomp(adjMat, 'Directed', false);
    exits.groupId = groupIds(:);
end