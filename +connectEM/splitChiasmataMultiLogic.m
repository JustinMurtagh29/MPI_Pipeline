function summary = splitChiasmataMultiLogic(summary)
    % summary = splitChiasmataMultiLogic(summary)
    %   This function decides which chiasmata can be safely split and which
    %   flight paths need to be applied to do so.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    chiCount = numel(summary.tracings);
    summary.split = false(chiCount, 1);
    
    for chiIdx = 1:numel(summary.tracings)
        tracings = summary.tracings{chiIdx};
        
        overlaps = reshape(tracings.overlaps, 1, []);
        overlaps = transpose(cell2mat(overlaps));
        nrExits = size(overlaps, 1);

        % partition chiasma
       [~, ~, lut] = buildPartition(overlaps);
       [partition, isValid, lut] = fixFourPartitions(overlaps, lut);
       
        if isequal(partition, [2; 2]) || ...
                (isValid && numel(partition) > 1)
            % do process
            % * all 2-2 partitions (valid and invalid)
            % * all valid partitions (except fully connected)
            splitChiasma = true;
            chiasmaIsSolved = false;
            executeTracings = selectTracings(overlaps);
        else
            splitChiasma = false;
            
            % valid fully connected chiasmata are solved
            chiasmaIsSolved = isValid && (numel(partition) == 1);
            executeTracings = false(nrExits, 1);
        end
        
        summary.split(chiIdx) = splitChiasma;
        summary.solved(chiIdx) = chiasmaIsSolved;
        summary.tracings{chiIdx}.execute = executeTracings;
        summary.tracings{chiIdx}.partition = partition;
        summary.tracings{chiIdx}.partitionIsValid = isValid;
        summary.tracings{chiIdx}.lut = lut;
    end
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

function [partition, isValid, lut] = buildPartition(overlaps)
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

function [partition, isValid, lut] = fixFourPartitions(overlaps, lut)
    % Email excerpt from amotta to L4team@mhlab.net (16.10.2017)
    % 
    % When we have flight paths A→B, B→A, C→D, D→A, this results in a 4-
    % partition. Assuming that there are no true branch points with 4
    % exits (which is not true, I've found one such bouton), the wrong
    % tracing can automatically be spotted (D→A) and corrected (D→C),
    % yielding a 2-2-partition. This error-correction is right in all
    % examples I've looked at.
    rawPartition = accumarray(lut, 1);
    fourPartitions = find(rawPartition == 4);
    
    for curPartition = reshape(fourPartitions, 1, [])
        curRows = find(lut == curPartition);
        curOverlaps = overlaps(curRows, :);
        
       [curUniOverlaps, ~, curUniCount] = ...
           unique(sort(curOverlaps, 2), 'rows');
        curUniOverlapsCount = accumarray(curUniCount, 1);
        curValidOverlap = curUniOverlaps(curUniOverlapsCount > 1, :);
        
        curErrorId = ismember(curOverlaps, curValidOverlap);
        curErrorId = find(xor(curErrorId(:, 1), curErrorId(:, 2)));
        
        % only proceed for valid 
        if numel(curErrorId) ~= 1; continue; end
        overlaps(curRows(curErrorId), 2) = 0;
    end
    
    % determine new partition
   [partition, isValid, lut] = buildPartition(overlaps);
end
