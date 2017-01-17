function enforceMyelinSegments(p, bbox, newPrefix)
    % enforceMyelinSegments(p, box)
    %   The current membrane detection / segmentation pipeline
    %   struggles with myelin sheaths and produces 
    %   neuron-myelin mergers. This function uses the results of
    %   myelin detection and patches up the previously generated
    %   classification result to reduce the rate of mergers.
    %
    % p
    %   Parameter structure from pipeline repo for dataset
    %
    % bbox
    %   Bounding Box where the fix is to be applied
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    % Adpoted as first child by:
    %   Manuel
 
    % make sure that the bounding box
    % is aligned with the KNOSSOS cube hierarchy
    assert(Util.checkBoundingBox(bbox));
    
    % Load myelin detection
    myelin = loadSegDataGlobal(p.myelin, bbox);
    border = bwdist(~myelin, 'euclidean');
    border = (border > 0) & (border < 2);

    % Load classification data
    class = loadClassData(p.class, bbox);
    
    % Build up the wall
    borderHeight = -1.7;
    class(border) = borderHeight;
    
    % write result
    saveClassData( ...
        Util.modifyStruct(param.class, 'prefix', newPrefix), bbox, class);
end

function [th5, th6] = buildMasks(data)
    % [plusMask, minusMask] = buildMasks(data)
    %   Very primitive function to locate large, black things.
    %   In practice, this primarily locates myelin sheaths and
    %   mitochondria.
    %
    % plusMask
    %   Grayscale image where positive values correspond to
    %   myelin sheaths and mitochondria.
    %
    % minusMask
    %   Grayscale image that highlights approx. 3 voxel thick
    %   borders around myelin and mitochondria.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % config
    thresh = 85;
    
    % 1 - threshold image
    th1 = (data > thresh);
    
    % 2 - remove shot noise
    th2 = imopen(th1, makeSphere(1));
    
    % 3 - remove membranes
    th3 = imdilate(th2, makeSphere(3));
    
    % 4 - remove small crap
    th4 = ~bwareaopen(~th3, 1E4);
    
    % 5 - restore whole myelin thickness
    % NOTE: myelin and mitos are now "true"
    th5 = ~th2 & (bwdist(~th4, 'euclidean') <= 5);
    th5 = bwmedian(th5, makeSphere(3));
    
    % 6 - find borders (thickness of two voxels)
    th6 = bwdist(~th5, 'euclidean');
    th6 = (th6 > 0) & (th6 <= 3);
end

function sphere = makeSphere(rad)
    % Based on code by Benedikt Staffler
    % From Util.getPointsInBall
    
    r2 = ceil(rad);
    [xx, yy, zz] = meshgrid(-r2:r2, -r2:r2, -r2:r2);
    sphere = sqrt(xx .^ 2 + yy .^ 2 + zz .^ 2) <= rad;
end

function out = bwmedian(im, strel)
    % out = bwmedian(im, strel)
    %   Fast n-dimensional median filter for BINARY images.
    %
    % im
    %   Input image: n-dimensional logical matrix
    %
    % strel
    %   Structuring element built with `strel` command
    
    out = convn(single(im), single(strel), 'same');
    out = out / sum(strel(:));
    out = (out > 0.5);
end

