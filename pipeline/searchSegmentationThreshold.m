function outFiles = searchSegmentationThreshold(p, bbox, tRange, doClass)
    % outFiles = searchSegmentationThreshold(p, bbox, tRange, doClass)
    %   Function that allows you to create multiple segmentation
    %   preview movies of different thresholds.
    %
    % p
    %   Parameter structure
    %
    % bbox
    %   Bounding box around region of interest
    %
    % tRange
    %   A vector of threshold values. For each of these values, a
    %   separate movie will be generated.
    %
    % doClass
    %   Flag which can be used to disable the classification step.
    %   If you're running this function multiple times on the same
    %   region of interest, you can set this flag to false.
    %   Default: true
    %
    % outFiles
    %   Cell array with the same shape as `tRange`. Each entry
    %   contains the path to the corresponding preview movie.
    %
    % Written by
    %   Jakob Straehle <jakob.straehle@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>

    if ~exist('doClass', 'var')
        doClass = true;
    end
    
    % prepare output
    tCount = numel(tRange);
    outFiles = cell(size(tRange));

    for curIdx = 1:tCount
        curT = tRange(curIdx);
        
        % change threshold and threshold function
        p.seg.threshold = curT;
        p.seg.func = @(x) watershedSeg_v1_cortex(x, {curT, 10});
        
        fprintf('Making movie %d of %d\n', curIdx, tCount);
        outFiles{curIdx} = makeSegmentationPreviewMovie(p, bbox, doClass);

        % prevent classification in subsequent runs
        doClass = false;
    end
end
