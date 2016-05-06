function [Wd,mask] = sparseKernel( W, d )
%SPARSEKERNEL Create d-regular sparse kernel from W.
% INPUT W: 4D array, where the fourth dimension corresponds to feature maps
%          and will not be transformed.
%       d: Vector containing the regularity for the first three dimensions
%          in W. (d = [2 2 2] means that every second row, colum, z-plane
%          will be zeros.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

fSize = size(W);
fdSize = fSize;
if length(fSize) < 3
    fSize(3) = 1;
end
fSize = fSize(1:3);
fsSize = fSize + (fSize - 1).*(d - 1);
fdSize(1:3) = fsSize;
Wd = zeros(fdSize,'like',W);
xg = false(fsSize(1),1);
yg = false(fsSize(2),1);
zg = false(fsSize(3),1);
xg(1:d(1):end) = true;
yg(1:d(2):end) = true;
zg(1:d(3):end) = true;
[x,y,z] = meshgrid(xg,yg,zg);
mask = x & y & z;
repSize = fdSize;
repSize(1:3) = 1;
Wd(repmat(mask,repSize)) = W;


end

