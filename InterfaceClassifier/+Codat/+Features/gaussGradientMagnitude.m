function [feat, border, numFeatures] = gaussGradientMagnitude(raw, sigma, truncate)
%GAUSSGRADIENTMAGNITUDE Multi-dimensional magnitude of gauss gradient.
% INPUT raw: N-dimensional input array.
%       sigma: Integer scalar or vector specifying the standard deviation of
%              the gaussian kernel for all or seperately for each dimension
%              of the input image.
%       truncate: (Optional) Double. Truncate filter at this many standard
%                 deviations.
%                 (Default: 3)
% OUTPUT feat: Array of the same size as raw containing the filter response.
%        border: Integer array specifying the total border of the filter
%                in each dimension.
%        numFeatures: Resulting feature space size per voxel.
%
% NOTE It is possible to calculate only the last two outputs by calling the
%      function without the raw input argument.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%parse inputs
if ~isempty(raw) && length(sigma) == 1
    sigma = repmat(sigma,[1, ndims(raw)]);
elseif ~isempty(raw) && length(sigma) ~= ndims(raw)
    error('Specify the standard deviation for all input dimensions.');
end
if ~exist('truncate','var') || isempty(truncate)
    truncate =  3;
end
border = 2.*ceil(sigma*truncate);
numFeatures = 1;
if isempty(raw)
    feat = [];
    return;
end

feat = zeros(size(raw),'like',raw);
for dim = 1:ndims(raw)
    order = zeros(ndims(raw),1);
    order(dim) = 1;
    feat = feat + Codat.Features.gaussianFilter(raw, sigma, truncate, order).^2;
end
feat = sqrt(feat);
end
