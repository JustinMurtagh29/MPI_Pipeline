function segmentation = watershedSeg_v2_cortex( aff, cell)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

tRange = cell{1};
vRange = cell{2};

segmentation = cell([length(tRange) length(vRange)]);

for t=1:length(tRange)
    fgm1 = aff < tRange(t);
    for v=1:length(vRange)
        fgm1 = bwareaopen(fgm1, vRange(v), 26);
        affImposed = imimposemin(aff, fgm1);
        segmentation{t,v} = watershed(affImposed, 26);
    end
end
end
