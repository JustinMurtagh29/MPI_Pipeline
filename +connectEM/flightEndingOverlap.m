function flightEndings = flightEndingOverlap( ...
        param, origAgglos, endings, flightNodes, flightAgglos, superAgglos)
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
    %   flightNodes
    %     Cell array with the nodes of each flight path.
    %
    %   flightAgglos
    %     Cell array with the IDs of super-agglomerates reached for each
    %     flight path. For each of the reached super-agglomerates up to one
    %     ending will be detected.
    %
    %   superAgglos
    %     Super-agglomerates to which the flight queries attach
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>

    assert(numel(flightNodes) == numel(flightAgglos));    
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
    flightEndings = cell(numel(flightNodes), 1);
    
    for curFlightIdx = 1:numel(flightEndings)
        fprintf('Processing flight #%d ...\n', curFlightIdx);
        
        % super-agglomerates reached by flight path
        curFlightAgglos = flightAgglos{curFlightIdx};
        curFlightAggloCount = numel(curFlightAgglos);
        
        % initialize output with zeros
        flightEndings{curFlightIdx} = zeros(curFlightAggloCount, 1);
        if ~curFlightAggloCount; continue; end;
        
        % convert flight path to nm
        curFlightNodes = double(flightNodes{curFlightIdx});
        curFlightNodes = bsxfun(@times, curFlightNodes, voxelSize);
        
        % for each super-agglomerate search ending
        for curFlightAggloIdx = 1:numel(curFlightAgglos)
            curFlightAgglo = curFlightAgglos(curFlightAggloIdx);
            
            % get endings in agglomerate of interest
            curEndings = aggloEndings(curFlightAgglo, :);
            curEndingCount = numel(curEndings{1});
            if ~curEndingCount; continue; end;
            
            % calculate minimum distance to endings
            curEndingDists = nan(curEndingCount, 1);
            for curEndingIdx = 1:curEndingCount
                curBorderPos = curEndings{2}{curEndingIdx};

                curDists = pdist2( ...
                    curFlightNodes, curBorderPos, 'euclidean');
                curEndingDists(curEndingIdx) = min(curDists(:));
            end
            
            % no ending reached if above distance threshold
           [curMinDist, curMinDistIdx] = min(curEndingDists);
            if curMinDist > maxDist; continue; end;
            
            % find ending with minimum distance
            curEndingId = curEndings{1}(curMinDistIdx);
            flightEndings{curFlightIdx}(curFlightAggloIdx) = curEndingId;
        end
    end
end

function aggloEndings = buildAggloEndings(voxelSize, endings, aggloOrigIds)
    %% group endings within original agglomerates
    origAggloCount = numel(endings.axonMask);
    origAggloEndings = cell(origAggloCount, 2);
    
    curEndingOff = 0;
    for curIdx = 1:numel(endings.axonIds)
        curAggloId = endings.axonIds(curIdx);
        curClusterIds = endings.borderClusters{curIdx};
        
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
