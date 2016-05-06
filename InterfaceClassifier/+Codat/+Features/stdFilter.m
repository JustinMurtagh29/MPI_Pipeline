function [feat, border, numFeatures] = stdFilter( raw, siz )
%MAXIMUMFILTER Calculates the local standard deviation of an input array.
% INPUT raw: n or 3 -dimensional input array.
%       siz: Integer vector or integer scalar specifying the size of the
%            filter sepeartely for each dimension or for all dimensions
%            respectively.
% OUTPUT feat: Array of the same size as raw containing the filter response.
%        border: Integer array specifying the total border of the filter
%                in each dimension.
%        numFeatures: Resulting feature space size per voxel.
%
% NOTE It is possible to calculate only the last two outputs by calling the
%      function without the raw input argument.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

border = siz - 1;
numFeatures = 1;
if isempty(raw)
    feat = [];
    return;
end
if length(siz) == 1
    siz = repmat(siz, ones(1, ndims(raw)));
end
feat = stdfilt(raw,ones(siz));
end

