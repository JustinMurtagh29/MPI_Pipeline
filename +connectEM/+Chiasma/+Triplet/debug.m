function skels = debug( ...
    origAxons, splitAxons, summary, flights, axonIds, chiasmaIds)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    skels = cell(numel(axonIds), 1);
    
    summaryAxonIds = cat(1, summary.axonId);
   [~, axonSummaryIds] = ismember(axonIds, summaryAxonIds);
    assert(all(axonSummaryIds));
    
    for curIdx = 1:numel(axonIds)
        curAxonId = axonIds(curIdx);
        
        curSummaryIdx = axonSummaryIds(curIdx);
        curSummary = summary(curSummaryIdx);
        
        % select chiasmata
        curChiasmaIds = chiasmaIds{curIdx};
        
        if isempty(curChiasmaIds)
            curChiasmaIds = 1:numel(curSummary.tracings);
        end
        
        curChiasmaIds = reshape(curChiasmaIds, [], 1);

       [curTaskIds, curOverlaps] = cellfun(@(t) ...
            deal(t.taskIds, cell2mat(transpose(t.overlaps))), ...
            curSummary.tracings(curChiasmaIds), ...
            'UniformOutput', false);
        curTaskIds = cat(1, curTaskIds{:});
        curOverlaps = cell2mat(curOverlaps')';
        
       [~, curFlightIds] = ismember( ...
           curTaskIds, flights.filenamesShort);
        assert(all(curFlightIds));
        
        curSkel = skeleton();
        
        % original axon
        curOrigAxon = origAxons(curAxonId);
        curSkel = curSkel.addTree( ...
            'Original', curOrigAxon.nodes(:, 1:3), curOrigAxon.edges);
        
        % split axons
        curSplitAxons = splitAxons{curAxonId};
        for curSplitIdx = 1:numel(curSplitAxons)
            curSplitAxon = curSplitAxons(curSplitIdx);
            
            curSkel = curSkel.addTree( ...
                sprintf('Split %d', curSplitIdx), ...
                curSplitAxon.nodes(:, 1:3), curSplitAxon.edges);
        end
        
        % flights
        for curFlightIdx = 1:numel(curFlightIds)
            curFlightId = curFlightIds(curFlightIdx);
            
            curTaskId = flights.filenamesShort{curFlightId};
            curFlightOverlaps = curOverlaps(curFlightIdx, :);
            
            curName = sprintf( ...
                'Flight %s (%d → %d)', ...
                curTaskId, curFlightOverlaps);
            
            curFlightNodes = flights.nodes{curFlightId};
            curSkel = curSkel.addTree(curName, flip(curFlightNodes, 1));
        end
        
        skels{curIdx} = curSkel;
    end
end