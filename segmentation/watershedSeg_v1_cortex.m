function segmentation = watershedSeg_v1_cortex( aff, cell )

    hRange = cell{1};                             
    vRange = cell{2};

    segmentation = cell([length(hRange) length(vRange)]);

    for h=1:length(hRange)
        affHmin= imhmin(aff, hRange(h), 26);  % suppress minima below hRange (def 0.25)
        bw1 = imregionalmin(affHmin, 26);   % mask with ones at regional minima
        clear affHmin;
        for v=1:length(vRange)
            bw1 = bwareaopen(bw1, vRange(v), 26);  % mask now only has objects bigger than vRange pixels (def 10)
            affImposed = imimposemin(aff, bw1);  % changes aff to remove minima that were smaller than vRange pixels
            segmentation{h,v} = watershed(affImposed, 26);  % watershed
            segmentation{h,v}(aff == -2) = 0;   % remove mask outliers (for mirrorPad)
        end
    end

end

