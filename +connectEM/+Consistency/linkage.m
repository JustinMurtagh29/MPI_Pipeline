function out = linkage(area)
    % rows are axons
    % columns are dendrites
    numAxons = size(area, 1);
    out = nan(numAxons - 1, 3);
    outRow = 0;
    
    maxClusterId = numAxons;
    clusterIds = 1:numAxons;
    clusters = num2cell(1:numAxons);
    
    lastGcRun = 0;
    gcDuration = 0;
    
    % WARNING(amotta): For performance reasons, the `w` and `m` matrices
    % use the following convention: row are dendrites, columns are axons!
    area = area - mean(area(:), 'omitnan');
    w = fillmissing(transpose(area), 'constant', 0);
    m = 1 - isnan(transpose(area));
    
    corr = transpose(w) * w;
    norm = transpose(m) * m;
    corr = corr ./ norm;
    
    corr(~norm) = nan;
    corr(1:(numAxons + 1):end) = nan;
    
    while true
       [maxCorr, maxIdx] = max(corr(:));
        if isnan(maxCorr); break; end
        
       [idxB, idxA] = ind2sub(size(corr), maxIdx);
        assert(idxA < idxB);
        
        outRow = outRow + 1;
        out(outRow, 1:2) = clusterIds([idxA, idxB]);
        out(outRow, 3) = maxCorr;
        
        cluster = [clusters{idxA}; clusters{idxB}];
        
        curW = area(cluster, :);
        curW = mean(curW, 1, 'omitnan');
        curM = 1 - isnan(curW);
        curW(~curM) = 0;
        
        % Update weights and mask
        w(:, idxA) = curW;
        w(:, idxB) = nan;
        m(:, idxA) = curM;
        m(:, idxB) = 0;
        
        % Update correlation
        curCorr = curW * w;
        curNorm = curM * m;
        curCorr = curCorr ./ curNorm;
        
        curCorr(~curNorm) = nan;
        curCorr(idxA) = nan;
        
        corr(idxA, :) = curCorr;
        corr(idxB, :) = nan;
        corr(:, idxA) = curCorr;
        corr(:, idxB) = nan;
        
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
            Util.log('GC pause');
            gcDuration = tic();
            
            % Garbage collection
            mask = ~isnan(clusterIds);
            clusterIds = clusterIds(mask);
            clusters = clusters(mask);
            w = w(:, mask);
            m = m(:, mask);
            corr = corr(mask, :);
            corr = corr(:, mask);
            
            gcDuration = toc(gcDuration);
            lastGcRun = tic();
        end
    end
    
    out = out(1:outRow, :);
end
