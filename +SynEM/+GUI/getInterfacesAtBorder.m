function [idx, touchesBorder] = getInterfacesAtBorder( data, border )
%GETINTERFACESATBORDER Get the indices of the interfaces at the border of a
% training cube.
% INPUT data: struct
%           Struct containing a training data cube.
%           (see calculateInterfacesFromSeg)
%       border: [3x1] int
%           Centroid distance from the cube boundary that is considered as
%           the border. Only interface surfaces with a voxel touching the
%           border and with a centroid within the border distance from the
%           cube border are considered border interfaces.
% OUTPUT idx: [Nx1] logical
%           Logical indices for the border in the border region.
%        touchesBorder: [Nx1] logical
%           Logical indices of those borders that have a voxel touching the
%           cube border.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

cen = cell2mat(data.centroidList);

touchesBorder = false(size(cen, 1), 1);

sz = size(data.segments);
border = border(:)';
for i = 1:length(data.interfaceSurfaceList)
    X = Util.indToSubMat(size(data.segments), ...
        data.interfaceSurfaceList{i}(:));
    touchesBorder(i) = any(X(:) == 1) | any(any(bsxfun(@eq, X, sz)));
end


idx = any(bsxfun(@le, cen, border), 2) | ...
      any(bsxfun(@ge, cen, sz - border), 2);

end

