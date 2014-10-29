function segmentation = watershedSeg_v3_paper( aff, hRange, vRange )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

segmentation = cell([length(hRange) length(vRange)]);
affHmin = cell(size(aff));
for h=1:length(hRange)
    for dir=1:length(aff)
        affHmin{dir} = imhmin(aff{dir}, hRange(h), 26);
    end
    affHminAll = (affHmin{1} + affHmin{2} + affHmin{3}) / 3;
    fgm1 = imregionalmin(affHminAll, 26);
    bgm1 = sobelEdge3D( aff{1}, aff{2}, aff{3}, 3, 3, [], 1);
    for v=1:length(vRange)
        fgm1 = bwareaopen(fgm1, vRange(v), 26);
        segmentation{h,v} = watershed_3D(affHmin{1}, affHmin{2}, affHmin{3}, fgm1 | bgm1);
    end
end

end

                                                     