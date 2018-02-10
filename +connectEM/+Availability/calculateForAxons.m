function targetClassAvail = ...
        calculateForAxons(param, blockData, saveDists, axonIds)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    voxelSize = param.raw.voxelSize;
    blockSize = blockData.info.param.blockSize;
    blockSizeNm = blockSize .* voxelSize;
    
    % linearize matrix
    targetClassAreas = reshape( ...
        blockData.targetClassAreas, ...
        numel(blockData.targetClasses), []);
    
    %% calculate (physical) coordinates of blocks
    blockCount = size(blockData.targetClassAreas);
    blockCount = blockCount(2:end);
    
    blockCenters = cell(1, 3);
   [blockCenters{:}] = ind2sub( ...
        blockCount, (1:size(targetClassAreas, 2))');
    blockCenters = cell2mat(blockCenters);

    % to physical unit
    blockCenters = blockCenters .* blockSizeNm;

    %% prepare output
    targetClassAvail = nan( ...
        size(targetClassAreas, 1), numel(saveDists));
    
    for curAxonIdx = 1:numel(axonIds)
        curAxonId = axonIds(curAxonIdx);
        
        % calculate physical coordinates for blocks of axon
        curAxonCenters = blockData.axonBlocks{curAxonId};
        curAxonCenters = curAxonCenters .* blockSizeNm;
        
        % calculate shortest distance to all blocks
        % NOTE(amotta): This is slow as hell...
        curDists = inf(1, size(targetClassAreas, 2));
        for curCenterIdx = 1:size(curAxonCenters, 1)
            curDists = min(curDists, pdist2( ...
                curAxonCenters(curCenterIdx, :), ...
                blockCenters, 'squaredeuclidean'));
        end

        curDists = sqrt(curDists);
        
        % accumulate target surface areas over distance
       [curDists, curSortIds] = sort(curDists, 'ascend');
        curCumTargetClassAreas = cumsum( ...
            targetClassAreas(:, curSortIds), 2);
        
        % save values for select distanecs
        curSaveIds = arrayfun(@(d) ...
            find(curDists < d, 1, 'last'), saveDists);
        targetClassAvail(curAxonId, :) = ...
            curCumTargetClassAreas(:, curSaveIds);
    end
end