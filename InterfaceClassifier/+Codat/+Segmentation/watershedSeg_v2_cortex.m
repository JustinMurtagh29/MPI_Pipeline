function segmentation = watershedSeg_v2_cortex( bm, tRange, vRange)
%WATERSHEDSEG_V2_CORTEX Segmentation procedure v2.
% Watershed segmentation on the boundary map with imposed minima based on
% a thresholded boundary (tRange) with small objects removed (vRange)
% INPUT bm: 3d array of single/dobule containing the boundary map.
%           Boundaries should have value 1 and intracellular space value 0.
%       tRange: [Nx1] array of double specifying the threshold on the
%           boundary map.
%       vRange: [Nx1] array of double specifying the size of object in
%           pixels which are removed before minima are imposed.
% OUTPUT segmentation: Cell array of size length(hRange) x length(vRange).
%           Each cell contains the segmentation for the corresponding
%           parameters.
% Author: Manuel Berning <manuel.berning@brain.mpg.de>
% Modified by Benedikt Staffler <benedikt.staffler@brain.mpg.de>

segmentation = cell([length(tRange) length(vRange)]);

for t=1:length(tRange)
    fgm1 = bm < tRange(t);
    for v=1:length(vRange)
        fgm1 = bwareaopen(fgm1, vRange(v), 26);
        affImposed = imimposemin(bm, fgm1);
        segmentation{t,v} = watershed(affImposed, 26);
    end
end
end
