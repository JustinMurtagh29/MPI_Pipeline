function exits = selectExits(chiasmaT, overlaps, maxNrQueries)
    % exits = selectExits(chiasmaT, overlaps, maxNrQueries)
    %   Selects which exits to query next based on detected chiasmata
    %   (`chiasmaT`) and partial answers (`overlaps`).
    %
    %   If `overlaps` is an empty array, all exits are considered for
    %   querying. This is also the default value for `overlaps`.
    %
    %   For each chiasma at most `maxNrQueries` are selected. By default,
    %   `maxNrQueries` is set to one. Set `maxNrQueries` to `inf` if you
    %   want to query all exits.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    %% default values
    if ~exist('overlaps', 'var')
        % by default, consider all exits
        overlaps = [];
    end
    
    if ~exist('maxNrQueries', 'var') || isempty(maxNrQueries)
        % set default value for `maxNrQueries`
        maxNrQueries = 1;
    end
    
    %% find exits to query
    chiasmaCount = size(chiasmaT, 1);
    chiasmaT.exitIds  = cell(chiasmaCount, 1);
    
    for curIdx = 1:size(chiasmaT, 1)
        curAggloId = chiasmaT.aggloId(curIdx);
        curChiasmaId = chiasmaT.chiasmaId(curIdx);
        
        curNrExits = chiasmaT.nrExits(curIdx);
        curOpenExits = 1:curNrExits;
        
        if ~isempty(overlaps)
            % remove queried or reached exits
            curOverlaps = overlaps{curAggloId}{curChiasmaId};
            curOpenExits = setdiff(curOpenExits, curOverlaps);
        end
        
        % select an open exits
        curNrQueries = min(numel(curOpenExits), maxNrQueries);
        curExits = curOpenExits(1:curNrQueries);
        
        chiasmaT.exitIds{curIdx} = curExits(:);
    end
    
    %% build output
    exitCounts = cellfun(@numel, chiasmaT.exitIds);
    
    exits = table;
    exits.aggloId = repelem(chiasmaT.aggloId, exitCounts);
    exits.chiasmaId = repelem(chiasmaT.chiasmaId, exitCounts);
    exits.exitId = cell2mat(chiasmaT.exitIds);
end