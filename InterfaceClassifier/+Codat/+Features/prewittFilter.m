function [feat, border, numFeatures] = prewittFilter(raw)
%PREWITTFILTER Calculate a multidimensional prewitt filter.
% INPUT raw: N-dimensional input array.
% OUTPUT feat: Array of the same size as raw containing the filter response.
%        border: Integer array specifying the total border of the filter
%                in each dimension.
%        numFeatures: Resulting feature space size per voxel.
%
% NOTE It is possible to calculate only the last two outputs by calling the
%      function without the raw input argument.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

border = 2;
numFeatures = 1;
if isempty(raw)
    feat = [];
    return
end
feat = zeros(size(raw), 'like', raw);
h1 = [1, 0, -1];
h2 = ones(1,3);
nD = ndims(raw);
for dim = 1:nD
    sz = ones(1, nD);
    sz(dim) = 3;
    h1 = reshape(h1,sz);
    tmp =  convn(raw,h1,'same');
    for rdim = setdiff(1:nD,dim)
        sz = ones(1, nD);
        sz(rdim) = 3;
        h2 = reshape(h2,sz);
        tmp = convn(tmp,h2,'same');
    end
    feat = feat + abs(tmp);
end
feat = feat./nD;
end
