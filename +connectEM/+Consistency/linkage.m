function [links, initDist] = linkage(rawArea)
    % rows are axons
    % columns are dendrites
    numAxons = size(rawArea, 1);
    links = nan(numAxons - 1, 3);
    outRow = 0;
    
    maxClusterId = numAxons;
    clusterIds = 1:numAxons;
    clusters = num2cell(1:numAxons);
    
    lastGcRun = 0;
    gcDuration = 0;
    
    area = rawArea;
    sq = @(a) a .* a;
    
    dist = squareform(pdist(area, ...
        @(a, b) mean(sq(a - b), 2, 'omitnan')));
    dist(1:(numAxons + 1):end) = nan;
    if nargout >= 2; initDist = dist; end
    
    while true
       [minDist, minIdx] = min(dist(:));
        if isnan(minDist); break; end
        
       [idxB, idxA] = ind2sub(size(dist), minIdx);
        assert(idxA < idxB);
        
        outRow = outRow + 1;
        links(outRow, 1:2) = clusterIds([idxA, idxB]);
        links(outRow, 3) = minDist;
        
        cluster = [clusters{idxA}; clusters{idxB}];
        clusterArea = mean(rawArea(cluster, :), 1, 'omitnan');
        clusterDist = mean(sq(area - clusterArea), 2, 'omitnan');
        
        % Update state
        area(idxA, :) = clusterArea;
        area(idxB, :) = nan;
        
        dist(idxA, :) = clusterDist;
        dist(:, idxA) = clusterDist;
        dist(idxA, idxA) = nan;
        dist(idxB, :) = nan;
        dist(:, idxB) = nan;
        
        clusters{idxA} = cluster;
        clusters{idxB} = zeros(0, 1);
        
        maxClusterId = maxClusterId + 1;
        clusterIds(idxA) = maxClusterId;
        clusterIds(idxB) = nan;
        
        % NOTE(amotta): For large axon counts the speed of the loop is
        % limited by `max` having to scan the entire `corr` matrix. Hence,
        % we should reduce the size of the correlation matrix (of the whole
        % problem, actually). But due to the overhead of creating the new,
        % smaller matrix, this shouldn't be done too often.
        if ~lastGcRun || toc(lastGcRun) > 10 * gcDuration
            gcDuration = tic();
            
            % Garbage collection
            mask = ~isnan(clusterIds);
            clusterIds = clusterIds(mask);
            clusters = clusters(mask);
            area = area(mask, :);
            dist = dist(mask, mask);
            
            gcDuration = toc(gcDuration);
            lastGcRun = tic();
        end
    end
    
    links = links(1:outRow, :);
    links(:, 1:2) = sort(links(:, 1:2), 2);
    
    % Add links that are missing due to undefined distance
    missingIds = setdiff(1:links(end, 2), links(:, 1:2));
    missingIds = reshape(missingIds, [], 1);
    
    missingLinks = reshape(links(end, 2) + (1:numel(missingIds)), [], 1);
    missingLinks = [missingIds, missingLinks, nan(size(missingIds))];
    links = [links; missingLinks];
end
