function [ bm ] = morphRecon( bm, r, range, scale )
%MORPHRECON Morphological reconstruction of boundary maps.
% INPUT bm: 3d array containing the boundary map.
%       r: Double radius for morphological operations.
%       range: [1x2] array of double. The first entry contains the value
%           corresponding to the boundary class and the second entry the
%           label for the intracellular class. bm values will be restricted
%           to the interval defined by range.
%       scale: (Optional) [1x3] array of double specifying the voxel
%           scaling which is used to determine the structuring element for
%           morphological reconstruction operations.
%           (Default: [1, 1, 1]);
% OUTPUT bm: Processed boundary maps with values in [0, 1], where 1
%           corresponds to the boundary class.
% Author: Manuel Berning <manuel.berning@brain.mpg.de>
% Modified by: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('scale','var') || isempty(scale)
    scale = [1, 1, 1];
end

if range(1) > range(2)
    bm = -bm;
    range = - range;
end

bm(bm > range(2)) = range(2);
bm(bm < range(1)) = range(1);
bm = (bm - range(1))./(range(2) - range(1));

if r ~= 0
    [x,y,z] = meshgrid(-r:r,-r:r,-r:r);
    se = (scale(1).*x/r).^2 + (scale(2).*y/r).^2 + (scale(3).*z/r).^2 <= 1;
    % Opening by reconstruction
    affEroded = imerode(bm, se);
    affRecon = imreconstruct(affEroded, bm);
    % Closing by reconstruction
    affReconDilated = imdilate(affRecon, se);
    bm = imreconstruct(imcomplement(affReconDilated), imcomplement(affRecon));
else
    bm = 1 - bm;
end


end

