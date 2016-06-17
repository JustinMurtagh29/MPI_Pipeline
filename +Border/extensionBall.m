function [idsOff, ballSize] = ...
        extensionBall(voxelSize, box, radiusInNm)
    % [idsOff, ballSize] = ...
    %     extensionBall(voxelSize, tileSize, radiusInNm)
    %   Computes the relative linear indices of all voxels
    %   contained in the ball with radius 'radiusInNm'.
    %
    % voxelSize
    %   Vector with the voxel size in nano-metres.
    %
    % box
    %    Bounding box in which the produced linear indices
    %    will be valid.
    %
    % radiuInNm
    %   Ball radius in nano-metres.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % compute radius in voxels
    boxSize = 1 + box(:, 2)' - box(:, 1)';
    radiusInVox = ceil(radiusInNm ./ voxelSize);
    
    % ATTENTION
    %   MATLAB f**ked up meshgrid
    [yMat, xMat, zMat] = meshgrid( ...
        voxelSize(2) .* ((-radiusInVox(2)):radiusInVox(2)), ...
        voxelSize(1) .* ((-radiusInVox(1)):radiusInVox(1)), ...
        voxelSize(3) .* ((-radiusInVox(3)):radiusInVox(3)));
    
    % build mask
    distMat = xMat .^ 2 + yMat .^ 2 + zMat .^ 2;
    mask = (distMat <= radiusInNm ^ 2);
    
    % build indices
    linIds = find(mask);
    ballSize = size(mask);
    
    [idsX, idsY, idsZ] = ...
        ind2sub(ballSize, linIds);
    
    % center indices on zero
    idsX = idsX - (radiusInVox(1) + 1);
    idsY = idsY - (radiusInVox(2) + 1);
    idsZ = idsZ - (radiusInVox(3) + 1);
    
    % build output
    idsOff = signedSubToInd( ...
        boxSize, [idsX, idsY, idsZ]);
    idsOff = int32(idsOff);
end

function ids = signedSubToInd(siz, subs)
    % ids = signedSubToInd(siz, subs)
    %   Converts signed, zero-based subscripts into
    %   linear indices. Yay!
    
    % row vector
    siz = siz(:)';
    
    % compute factor for each dimension
    dimFactor = [1, siz(1:(end - 1))];
    dimFactor = cumprod(dimFactor);
    
    % build output
    ids = bsxfun(@times, subs, dimFactor);
    ids = sum(ids, 2);
end