function chiasma = selectExits(agglos, chiasmata, overlaps)
    % chiasma = selectExits(agglos, chiasmata, overlaps)
    %   Selects which exits to query next based on detected chiasmata
    %   (`chiasmata`) and partial answers (`overlaps`). Chiasmata marked as
    %   solved (`agglos.solvedChiasma`) are ignored.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    assert(numel(agglos) == numel(chiasmata));
    assert(numel(agglos) == numel(overlaps));
    
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
    chiasma.exitId  = zeros(chiasmaCount, 1);
    
    for curIdx = 1:size(chiasma, 1)
        curAggloId = chiasma.aggloId(curIdx);
        curChiasmaId = chiasma.chiasmaId(curIdx);
        
        curChiasmata = chiasmata{curAggloId};
        curNrExits = numel(curChiasmata.queryIdx{curChiasmaId});
        curOverlaps = overlaps{curAggloId}{curChiasmaId};
        
        % select an open ending
        curOpenExits = setdiff(1:curNrExits, curOverlaps);
        curExit = max(cat(1, 0, min(curOpenExits)));
        
        chiasma.exitId(curIdx) = curExit;
    end
end