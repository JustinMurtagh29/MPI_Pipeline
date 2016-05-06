function [feat, border, numFeatures] = differenceOfGaussians(raw, sigma, k, truncate)
%DIFFERENCEOFGAUSSIANS Multidimensional difference of gaussians filter.
% INPUT raw: N-dimensional input array.
%       sigma: Integer scalar or vector specifying the standard deviation of
%           the fist gaussian kernel for all or seperately for each
%           dimension of the input image. For the standard deviation of the
%           second gaussian see input 'k'.
%       k: (Optional) The standard deviation of the second gaussian is the
%          standard deviation of the first gaussian multiplied by the
%          scalar k.
%          (Default: 1.5)
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
if ~exist('k','var') || isempty(k)
    k =  1.5;
end
if ~exist('truncate','var') || isempty(truncate)
    truncate = 3;
end
border = 2*ceil(k*sigma*truncate);
numFeatures = 1;
if isempty(raw)
    feat = [];
    return;
end

I1 = Codat.Features.gaussianFilter(raw, sigma, truncate);
I2 = Codat.Features.gaussianFilter(raw, k*sigma, truncate);

feat = I1 - I2;
end
