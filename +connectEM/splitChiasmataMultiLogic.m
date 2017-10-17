function [splitChiasma, executeFlights] = splitChiasmataMultiLogic(summary)
    % [splitChiasma, executeFlights] = splitChiasmataMultiLogic(summary)
    %   This function decides which chiasmata can be safely split and which
    %   flight paths need to be applied to do so.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
   [splitChiasma, executeFlights] = cellfun( ...
       @forChiasma, summary.tracings, 'UniformOutput', false);
    splitChiasma = cell2mat(splitChiasma);
end

function [splitChiasma, executeTracings] = forChiasma(tracings)
    overlaps = reshape(tracings.overlaps, 1, []);
    overlaps = transpose(cell2mat(overlaps));
    nrExits = size(overlaps, 1);

    splitChiasma = false;
    executeTracings = false(nrExits, 1);
    
    % partition chiasma
   [partition, isValid] = buildPartition(overlaps);
   
    % do not cut out sphere if there's nothing we can do
    if isequal(partition, nrExits); return; end
    
    % do not process chiasmata with invalid solutions (for now)
    if ~isValid; return; end
    
    % select flight paths
    splitChiasma = true;
    executeTracings = selectTracings(overlaps);
end

function execute = selectTracings(overlaps)
    overlaps = sort(overlaps, 2);
    nrExits = size(overlaps, 1);
    
    adjMat = overlaps;
    adjMat(~all(adjMat, 2), :) = [];
    
    % build minimal spanning tree per component
    adjMat = sparse(adjMat(:, 2), adjMat(:, 1), 1, nrExits, nrExits);
   	adjMat = graphminspantree(adjMat, 'Method', 'Kruskal');
    
   [edges(:, 2), edges(:, 1)] = find(adjMat);
   
    % find flights which form the minimal spanning tree
   [~, tracingIds] = ismember(edges, overlaps, 'rows');
   
    % build mask
    execute = false(nrExits, 1);
    execute(tracingIds) = true;
end

function [partition, isValid] = buildPartition(overlaps)
    % build effective edges
    edges = sort(overlaps, 2);
    edges(~all(edges, 2), :) = [];
    edges = unique(edges, 'rows');
    edges = reshape(edges, [], 2);
    
    % find parition
    nrExits = size(overlaps, 1);
    adjMat = sparse(edges(:, 2), edges(:, 1), true, nrExits, nrExits);
   [~, lut] = graphconncomp(adjMat, 'Directed', false);
    lut = reshape(lut, [], 1);
   
    % build partition
    partition = accumarray(lut, 1);
    partition = sort(partition, 'descend');
    
    % For each component with at least two elements, there exists an edge.
    % Hence, there is no excuse of the flight paths involved in these
    % components to be dangling.
    isValid = ~any( ...
        overlaps(:, 2) == 0 & ismember( ...
        overlaps(:, 1), overlaps(:, 2)));
end