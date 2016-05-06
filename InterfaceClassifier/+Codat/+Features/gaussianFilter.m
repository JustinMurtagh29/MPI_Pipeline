function [feat, border, numFeatures] = gaussianFilter(raw, sigma, truncate, order)
%GAUSSIANFILTER Multidimensional gaussian filter.
% INPUT raw: N-dimensional input array.
%       sigma: Integer scalar or vector specifying the standard deviation of
%              the gaussian kernel for all or seperately for each dimension
%              of the input image. A value of 0 corresponds to no filtering
%              along this dimension.
%       truncate: (Optional) Double. Truncate filter at this many standard
%                 deviations.
%                 (Default: 3)
%       order: (Optional) Integer vector specifying the degree of the
%           derivative of the gaussian filter in each dimension. An order
%           of 0 corresponds to the convolution with a gaussian kernel and
%           orders of 1 or 2 correspond to convolutions with the first or
%           second derivative of the gaussian kernel. Higher orders are not
%           implemented.
%           (Default: 0)
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
if ~exist('order','var') || isempty(order)
    order =  zeros(length(sigma), 1);
end
border = 2.*ceil(sigma*truncate);
numFeatures = 1;
if isempty(raw)
    feat = [];
    return;
end

feat = raw;
for dim = 1:ndims(raw)
    if sigma(dim) == 0
        continue;
    end
    coords = ceil(sigma(dim)*truncate);
    sz = ones(1,length(sigma));
    sz(dim) = 2*coords + 1;
    x = reshape(-coords:coords, sz);
    h = exp(-x.^2./2./sigma(dim).^2);
    h = h./sum(h(:));
    switch order(dim)
        case 0
            %already calculated
        case 1
            h = -x./sigma(dim)^2.*h;
        case 2
            h = (x.^2./sigma(dim)^2 - 1)./sigma(dim).^2.*h;
        otherwise
            error('An order of ''%d'' is not implemented.', order(dim));
    end
    feat = convn(feat,h,'same');
end
end
