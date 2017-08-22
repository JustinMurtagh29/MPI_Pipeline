function flightEndings = flightEndingOverlap( ...
        param, origAgglos, endings, flights, flightResults, superAgglos)
    % flightEndingOverlap
    % 
    % Inputs
    %   param
    %     Parameter structure of pipeline run
    %
    %   origAgglos
    %     Original agglomerates on which the ending detection ran
    %
    %   endings
    %     Endings detected on the original agglomerates. For all but the
    %     first round of flight queries you might want to pass in the
    %     subset of remaining open endings.
    %
    %   flights
    %     Focused flights. Also known as the "ff" structure.
    %
    %   flightResults
    %     Results of the focused flight analysis. It is expected that all
    %     flights passed into this function have a valid end agglomerate.
    %
    %   superAgglos
    %     Super-agglomerates to which the flight queries attach
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    maxSegId = Seg.Global.getMaxSegId(param);
    voxelSize = param.raw.voxelSize;
    maxDist = 300; % in nm
    
    % build look-up table
    % maps segment ID → original agglomerate ID
    origAggloLUT = buildLUT(maxSegId, origAgglos);
    
    % convert super-agglomerates to regular agglomerates
    agglos = arrayfun(@(sa) {sa.nodes(:, 4)}, superAgglos);
    agglos = cellfun(@(ids) {unique(ids(~isnan(ids)))}, agglos);
    
    % group endings
    aggloOrigIds = cellfun( ...
        @(ids) {setdiff(origAggloLUT(ids), 0)}, agglos);
    aggloEndings = buildAggloEndings(voxelSize, endings, aggloOrigIds);
    
    %% do the magic
    flightCount = numel(flights.nodes);
    flightEndings = zeros(flightCount, 1);
    
    for curFlightIdx = 1:flightCount
        fprintf('Processing flight #%d ...\n', curFlightIdx);
        curEndAggloId = flightResults.endAgglo{curFlightIdx};
        
        % NOTE(amotta): This should not be needed because the fight paths
        % are expected to be filtered before calling this function.
        if isempty(curEndAggloId); continue; end;
        
        if numel(curEndAggloId) > 1
            % HACK(amotta): This should be handled in the query analysis
            warning('Discarding endings for flight #%d', curFlightIdx)
            curEndAggloId = max(curEndAggloId);
        end
        
        curEndings = aggloEndings(curEndAggloId, :);
        curEndingCount = numel(curEndings{1});
        if ~curEndingCount; continue; end;
        
        % convert flight path to nm space
        curFlightNodes = flights.nodes{curFlightIdx};
        curFlightNodes = bsxfun(@times, curFlightNodes, voxelSize);
        
        curEndingDists = nan(curEndingCount, 1);
        for curEndingIdx = 1:curEndingCount
            curBorderPos = curEndings{2}{curEndingIdx};
            
            curDists = pdist2( ...
                curFlightNodes, curBorderPos, 'squaredeuclidean');
            curEndingDists(curEndingIdx) = sqrt(min(curDists(:)));
        end
        
       [curMinDist, curMinDistIdx] = min(curEndingDists);
        if curMinDist > maxDist; continue; end;
        
        curEndingId = curEndings{1}(curMinDistIdx);
        flightEndings(curFlightIdx) = curEndingId;
    end
end

function aggloEndings = buildAggloEndings(voxelSize, endings, aggloOrigIds)
    %% group endings within original agglomerates
    origAggloCount = numel(endings.idxCanidateFound);
    origAggloEndings = cell(origAggloCount, 2);
    
    curEndingOff = 0;
    for curIdx = 1:numel(endings.mapping)
        curAggloId = endings.mapping(curIdx);
        curClusterIds = endings.T{curIdx};
        
        % convert borders to nm space
        curBorderPos = double(endings.borderPositions{curIdx});
        curBorderPos = bsxfun(@times, curBorderPos, voxelSize);
        
        curEndingCount = max(curClusterIds);
        curEndingIds = curEndingOff + (1:curEndingCount);
        curEndingOff = curEndingOff + curEndingCount;
        
        % group borders
        curClusters = accumarray( ...
            curClusterIds, 1:numel(curClusterIds), ...
            [], @(rows) {curBorderPos(rows, :)});
        
        % save result
        origAggloEndings{curAggloId, 1} = curEndingIds(:);
        origAggloEndings{curAggloId, 2} = curClusters;
    end
    
    %% group for super agglomerates
    aggloCount = numel(aggloOrigIds);
    aggloEndings = cell(aggloCount, 2);
    
    for curIdx = 1:aggloCount
        curEndings = origAggloEndings(aggloOrigIds{curIdx}, :);
        
        aggloEndings{curIdx, 1} = cell2mat(curEndings(:, 1));
        aggloEndings{curIdx, 2} = cat(1, curEndings{:, 2});
    end
end

function lut = buildLUT(maxSegId, agglos)
    lut = zeros(maxSegId, 1);
    lut(cell2mat(agglos)) = repelem( ...
        1:numel(agglos), cellfun(@numel, agglos));
end