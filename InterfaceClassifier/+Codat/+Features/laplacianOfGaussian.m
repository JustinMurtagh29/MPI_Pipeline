function [feat, border, numFeatures] = laplacianOfGaussian(raw, sigma, truncate)
%LAPLACIANOFGAUSSIAN Laplace of gaussians filter.
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
    truncate = 3;
end
border = 2.*ceil(sigma*truncate);
numFeatures = 1;
if isempty(raw)
    feat = [];
    return;
end

feat = zeros(size(raw),'like',raw);
nD = ndims(raw);
for dim = 1:nD
    order = zeros(nD,1);
    order(dim) = 2;
    feat = feat + Codat.Features.gaussianFilter(raw, sigma, truncate, order);
end
end
