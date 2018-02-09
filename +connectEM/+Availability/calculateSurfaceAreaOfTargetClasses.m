function [classSurfAreas, classes] = ...
        calculateSurfaceAreaOfTargetClasses(conn, aggloSurfAreas)
    % Written by
    %   Kevin M. Boergens <kevin.boergens@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    %% build target classes
    classes = unique(conn.denMeta.targetClass);
    classes = sort(reordercats(classes));
    
   [~, classAggloIds] = ismember( ...
        conn.denMeta.targetClass, classes);
    classAggloIds = accumarray( ...
        classAggloIds, conn.denMeta.id, ...
        size(classes), @(ids) {ids(:)});
    
    % add extra category with all postsynaptic agglomerates
    allAggloIds = reshape(1:size(conn.dendrites), [], 1);
    classAggloIds = cat(1, {allAggloIds}, classAggloIds);
    classes = cat(1, 'All', classes);
    
    %% calculate per-class surface area
    classCount = numel(classes);
    blockCount = numel(aggloSurfAreas);
    
    classSurfAreas = nan(classCount, blockCount);
    for curBlockIdx = 1:blockCount
        curBlock = aggloSurfAreas{curBlockIdx};
        
        for curClassIdx = 1:classCount
            curSurfArea = classAggloIds{curClassIdx};
            curSurfArea = ismember(curBlock(:, 1:2), curSurfArea);
            curSurfArea = sum(curBlock(:, 4) .* sum(curSurfArea, 2));
            classSurfAreas(curClassIdx, curBlockIdx) = curSurfArea;
        end
    end
    
    classSurfAreas = reshape( ...
        classSurfAreas, [classCount, size(aggloSurfAreas)]);
end