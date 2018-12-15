function out = linkage(area)
    % rows are axons
    % columns are dendrites
    numAxons = size(area, 1);
    out = nan(numAxons - 1, 3);
    outRow = 1;
    
    maxClusterId = numAxons;
    clusterIds = 1:numAxons;
    clusters = num2cell(1:numAxons);
    
    area = area - mean(area(:), 'omitnan');
    w = fillmissing(area, 'constant', 0);
    m = 1 - isnan(area);
    
    corr = w * w';
    norm = m * m';
    corr = corr ./ norm;
    
    corr(~norm) = nan;
    corr(1:(numAxons + 1):end) = nan;
    
    while true
       [maxCorr, maxIdx] = max(corr(:));
        if isnan(maxCorr); break; end
        
       [idxB, idxA] = ind2sub(size(corr), maxIdx);
        assert(idxA < idxB);
        
        out(outRow, 1:2) = clusterIds([idxA, idxB]);
        out(outRow, 3) = maxCorr;
        outRow = outRow + 1;
        
        cluster = [clusters{idxA}; clusters{idxB}];
        
        curW = area(cluster, :);
        curW = mean(curW, 1, 'omitnan');
        curM = 1 - isnan(curW);
        curW(~curM) = 0;
        
        % Update weights and mask
        w(idxA, :) = curW;
        w(idxB, :) = nan;
        m(idxA, :) = curM;
        m(idxB, :) = 0;
        
        % Update correlation
        curCorr = w * curW(:);
        curNorm = m * curM(:);
        curCorr = curCorr ./ curNorm;
        
        curCorr(~curNorm) = nan;
        curCorr(idxA) = nan;
        
        corr(idxA, :) = curCorr;
        corr(idxB, :) = nan;
        corr(:, idxA) = curCorr;
        corr(:, idxB) = nan;
        
        clusters{idxA} = cluster;
        maxClusterId = maxClusterId + 1;
        clusterIds(idxA) = maxClusterId;
        clusterIds(idxB) = nan;
        
        % NOTE(amotta): For large axon counts, the speed of the loop is
        % limited by `max` having to scan the entire `corr` matrix. Hence,
        % we should reduce the size of the correlation matrix (of the whole
        % problem, actually). But due to the overhead of creating the new,
        % smaller matrix, this shouldn't be done too often.
        mask = isnan(clusterIds);
        if mean(mask) > 0.05
            clusterIds(mask) = [];
            clusters(mask) = [];
            w(mask, :) = [];
            m(mask, :) = [];
            corr(mask, :) = [];
            corr(:, mask) = [];
        end
    end
    
end
