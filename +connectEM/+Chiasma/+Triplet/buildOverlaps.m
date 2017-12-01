function overlaps = buildOverlaps(chiasmata, summary)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    overlaps = cellfun(@(a) ...
        arrayfun(@(~) ...
            zeros(3, 2), a.ccNodeIdx, ...
            'UniformOutput', false), ...
        chiasmata, 'UniformOutput', false);
    
    for curAxonIdx = 1:numel(summary)
        curSummary = summary(curAxonIdx);
        
        curAxonId = curSummary.axonId;
        curChiasmaIds = curSummary.chiasmaId;
        
        for curChiIdx = 1:numel(curChiasmaIds)
            curChiId = curChiasmaIds(curChiIdx);
            
            curOverlaps = curSummary.tracings{curChiIdx}.overlaps;
            curOverlaps = cell2mat(reshape(curOverlaps, 1, []));
            curOverlaps = transpose(curOverlaps);
            
            % complete overlaps
            curFull = zeros(3, 2);
            curFull(curOverlaps(:, 1), :) = curOverlaps;
            
            curOpenIds = find(~curFull(:, 1));
            curOpenIds = reshape(curOpenIds, 1, []);
            
            for curId = curOpenIds
                % check if there is an incoming flight
                curInId = find(curFull(:, 2) == curId, 1);
                if isempty(curInId); continue; end
                
                % if yes, use it to fake reverse flight
                curFull(curId, :) = cat(2, curId, curInId);
            end
            
            % check if there is a dangling flight
            curDangId = find(curFull(:, 1) & ~curFull(:, 2));
            assert(numel(curDangId) < 2);
            
            if isscalar(curDangId)
                % If there's a dangling flight path (e.g., originating from
                % exit A), we known / assume that B and C must be mutually
                % connected.
                curOtherIds = setdiff(1:3, curDangId);
                curFull(curOtherIds, 1) = curOtherIds;
                curFull(curOtherIds, 2) = fliplr(curOtherIds);
            end
            
            overlaps{curAxonId}{curChiId} = curFull;
        end
    end
end