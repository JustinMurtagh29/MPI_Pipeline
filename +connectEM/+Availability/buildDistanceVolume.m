function mask = buildDistanceVolume( ...
        param, distRefIso, distThreshUm, bbox)
	% Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    distThresh = 1E3 * distThreshUm;
    voxelSize = param.raw.voxelSize;
    
    boxCenter = mean(bbox, 2)';
    boxCenter = boxCenter .* voxelSize;
    
    boxSize = 1 + diff(bbox, 1, 2)';    
    boxRadius = boxSize .* voxelSize;
    boxRadius = sqrt(sum(boxRadius .^ 2)) / 2;
    
    distRefIso = distRefIso.vertices;
    distRefIso = distRefIso .* voxelSize;
    
    distRefDist = pdist2(distRefIso, boxCenter);
    if any((distRefDist + boxRadius) < distThresh)
        mask = true(boxSize);
        return;
    end
    
    mask = false(boxSize);
    distRefDist = (distRefDist - boxRadius) < distThresh;
    distRefIso = distRefIso(distRefDist, :);
    if isempty(distRefIso); return; end
    
    threshOne = distThresh;
    threshTwo = threshOne / sqrt(3);
    
    curDists = distRefIso;
    curMaskOne = false(size(curDists));
    curMaskTwo = false(size(curDists));

    for curX = 1:size(mask, 1)
        curDists(:, 1) = abs(distRefIso(:, 1) - ...
            (bbox(1, 1) + curX - 1) .* voxelSize(1));
        curMaskOne(:, 1) = curDists(:, 1) < threshOne;
        curMaskTwo(:, 1) = curDists(:, 1) < threshTwo;
        
        for curY = 1:size(mask, 2)
            curDists(:, 2) = abs(distRefIso(:, 2) ...
                - (bbox(2, 1) + curY - 1) .* voxelSize(2));
            curMaskOne(:, 2) = curMaskOne(:, 1) ...
                & (curDists(:, 2) < threshOne);
            curMaskTwo(:, 2) = curMaskTwo(:, 1) ...
                & (curDists(:, 2) < threshTwo);
            
            for curZ = 1:size(mask, 3)
                curDists(:, 3) = abs(distRefIso(:, 3) ...
                    - (bbox(3, 1) + curZ - 1) .* voxelSize(3));
                
                if any(curMaskTwo(:, 2) & (curDists(:, 3) < threshTwo))
                    mask(curX, curY, curZ) = true;
                    continue;
                end
                
                curMaskOne(:, 3) = curMaskOne(:, 2) ...
                    & (curDists(:, 3) < threshOne);
                if ~any(curMaskOne(:, 3)); continue; end
                
                curDist = curDists(curMaskOne(:, 3), :);
                curDist = sqrt(min(sum(curDist .* curDist, 2)));
                mask(curX, curY, curZ) = curDist < threshOne;
            end
        end
    end
end