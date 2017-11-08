function exits = selectExits(agglos, chiasmata, overlaps, maxNrQueries)
    % exits = selectExits(agglos, chiasmata, overlaps, maxNrQueries)
    %   Selects which exits to query next based on detected chiasmata
    %   (`chiasmata`) and partial answers (`overlaps`). Chiasmata marked as
    %   solved (`agglos.solvedChiasma`) are ignored.
    %
    %   For each chiasma at most `maxNrQueries` are selected. By default,
    %   `maxNrQueries` is set to one. Set `maxNrQueries` to `inf` if you
    %   want to query all exits.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    assert(numel(agglos) == numel(chiasmata));
    assert(numel(agglos) == numel(overlaps));
    
    if ~exist('maxNrQueries', 'var') || isempty(maxNrQueries)
        % set default value for `maxNrQueries`
        maxNrQueries = 1;
    end
    
    chiasma = table;
    chiasma.aggloId = repelem((1:numel(chiasmata))', ...
        cellfun(@(c) numel(c.ccCenterIdx), chiasmata));
    chiasma.chiasmaId = cell2mat(cellfun(@(c) ...
        (1:numel(c.ccCenterIdx))', chiasmata, 'UniformOutput', false));
    
    % get rid of solved chiasmata
    chiasma.isSolved = arrayfun(@(a, c) ...
        agglos(a).solvedChiasma(chiasmata{a}.ccCenterIdx(c)), ...
        chiasma.aggloId, chiasma.chiasmaId);
    chiasma(chiasma.isSolved, :) = [];
    chiasma.isSolved = [];
    
    %% find exits to query
    chiasmaCount = size(chiasma, 1);
    chiasma.exitIds  = cell(chiasmaCount, 1);
    
    for curIdx = 1:size(chiasma, 1)
        curAggloId = chiasma.aggloId(curIdx);
        curChiasmaId = chiasma.chiasmaId(curIdx);
        
        curChiasmata = chiasmata{curAggloId};
        curNrExits = numel(curChiasmata.queryIdx{curChiasmaId});
        curOverlaps = overlaps{curAggloId}{curChiasmaId};
        
        % select an open exits
        curOpenExits = setdiff(1:curNrExits, curOverlaps);
        curNrQueries = min(numel(curOpenExits), maxNrQueries);
        curExits = curOpenExits(1:curNrQueries);
        
        chiasma.exitIds{curIdx} = curExits(:);
    end
    
    %% build output
    exitCounts = cellfun(@numel, chiasma.exitIds);
    
    exits = table;
    exits.aggloId = repelem(chiasma.aggloId, exitCounts);
    exits.chiasmaId = repelem(chiasma.chiasmaId, exitCounts);
    exits.exitId = cell2mat(chiasma.exitIds);
    
end