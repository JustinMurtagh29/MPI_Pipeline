function distVol = buildDistanceVolume(param, distRefPoints, bbox)
    % distVol = buildDistanceVolume(param, distRefIso, bbox)
    %   This is extremely stupid!
    %
	% Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    voxelSize = param.raw.voxelSize;
    distRefPoints = distRefPoints .* voxelSize;
    boxSize = 1 + diff(bbox, 1, 2)';
    
    distVol = nan(boxSize);
    for curX = 1:size(distVol, 1)
        curDists(:, 1) = sq(distRefPoints(:, 1) - ...
            (bbox(1, 1) + curX - 1) .* voxelSize(1));
        
        for curY = 1:size(distVol, 2)
            curDists(:, 2) = ...
                curDists(:, 1) + sq(distRefPoints(:, 2) ...
                - (bbox(2, 1) + curY - 1) .* voxelSize(2));
            
            for curZ = 1:size(distVol, 3)
                curDists(:, 3) = ...
                    curDists(:, 2) + sq(distRefPoints(:, 3) ...
                    - (bbox(3, 1) + curZ - 1) .* voxelSize(3));
                distVol(curX, curY, curZ) = sqrt(min(curDists(:, 3)));
            end
        end
    end
end

function d = sq(d)
    d = d .* d;
end
