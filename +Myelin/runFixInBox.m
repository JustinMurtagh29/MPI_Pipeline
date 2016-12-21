function runFixInBox(param, newPrefix, box)
    % runFixInBox(param, newPrefix, box)
    %   The current membrane detection / segmentation pipeline
    %   struggles with myelin sheaths and produces a lot of
    %   neuron-myelin mergers. This function tries to identify
    %   myelin and patches up the previously generated classi-
    %   fication result to reduce the rate of mergers.
    %
    % NOTE
    %   * This is an ugly hack and was developped specifically
    %     for Jakob's mouse M1 L6 dataset (s40_L6M1_js_v10).
    %     Most constants are hard-coded for this very dataset.
    %
    %   * The algorithm for myelin detection is very primitive
    %     and produces false positives for mitochondria. This
    %     is not a problem as such, as the merge rate decreases.
    %     But this means that you *SHOULD NOT* use the output
    %     of this function to check whether a given axon is
    %     myelinated.
    %
    % param
    %   Parameter structure for current dataset
    %
    % newPrefix
    %   Prefix of the patched up classification results. The
    %   KNOSSOS cubes are written to the same root directory
    %   as the original classification.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % config
    pad = 10;

    % make sure that the bounding box
    % is aligned with the KNOSSOS cube hierarchy
    assert(Util.checkBoundingBox(box));
    
    % add padding
    disp('> Loading raw data...');
    boxBig = [box(:, 1) - pad, box(:, 2) + pad];
    raw = loadRawData(param.raw.root, param.raw.prefix, boxBig);
    
    % build masks around myelin (and mitos)
    disp('> Searching for myelin...');
    [plusMask, minusMask] = buildMasks(raw);
    clear raw;
    
    % remove padding
    plusMask = dropPadding(plusMask, pad);
    minusMask = dropPadding(minusMask, pad);
    
    disp('> Loading classification data...');
    class = loadClassData(param.class, box);
    
    % find a suitably small / big value
    smallBig = prctile(class(:), [5, 95]);
    % class(plusMask) = smallBig(end);
    class(minusMask) = smallBig(1);
    
    % write result
    disp('> Saving modified classification...');
    saveClassData(param.class.root, newPrefix, box, class);
end

function data = dropPadding(data, pad)
    data = data( ...
        (1 + pad):(end - pad), ...
        (1 + pad):(end - pad), ...
        (1 + pad):(end - pad));
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
