function segmentation = watershedSeg_v1_cortex( bm, hRange, vRange )
%WATERSHEDSEG_V1_CORTEX Segmentation procedure v1.
% Watershed segmentation on the boundary map with imposed minima based on
% regional minima with suppression (hRange) and removal of small objects
% (vRange).
% INPUT bm: 3d array of single/dobule containing the boundary map.
%           Boundaries should have value 1 and intracellular space value 0.
%       hRange: [Nx1] array of double containing the values used for minima
%           suppression.
%       vRange: [Nx1] array of double specifying the size of object in
%           pixels which are removed before minima are imposed.
% OUTPUT segmentation: Cell array of size length(hRange) x length(vRange).
%           Each cell contains the segmentation for the corresponding
%           parameters.
% Author: Manuel Berning <manuel.berning@brain.mpg.de>
% Modified by Benedikt Staffler <benedikt.staffler@brain.mpg.de>

segmentation = cell([length(hRange) length(vRange)]);

for h=1:length(hRange)
    affHmin = imhmin(bm, hRange(h), 26);
    bw1 = imregionalmin(affHmin, 26);
    clear affHmin;
    for v=1:length(vRange)
        bw1 = bwareaopen(bw1, vRange(v), 26);
        affImposed = imimposemin(bm, bw1);
        segmentation{h,v} = watershed(affImposed, 26);
    end
end

end