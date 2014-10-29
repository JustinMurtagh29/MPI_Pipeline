function segmentation = watershedSeg_v1_cortex( aff, cell )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

hRange = cell{1};                             
vRange = cell{2};

segmentation = cell([length(hRange) length(vRange)]);

for h=1:length(hRange)
    affHmin= imhmin(aff, hRange(h), 26);
    bw1 = imregionalmin(affHmin, 26);
    clear affHmin;
    for v=1:length(vRange)
        bw1 = bwareaopen(bw1, vRange(v), 26);
        affImposed = imimposemin(aff, bw1);
        segmentation{h,v} = watershed(affImposed, 26);
    end
end

end