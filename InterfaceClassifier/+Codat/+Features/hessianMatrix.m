function [H, border, numFeatures] = hessianMatrix(raw, sigma, truncate)
%HESSIANMATRIX Subsumed second derivatives of N-dimensional array.
% INPUT raw: N-dimensional input array.
%       sigma: Integer scalar or vector specifying the standard deviation of
%              the gaussian kernel for all or seperately for each dimension
%              of the input image.
%       truncate: (Optional) Double. Truncate filter at this many standard
%                 deviations.
%                 (Default: 3)
% OUTPUT H: Square cell array of size ndims(raw) in both dimension. Each
%           cell contains the second derivative with respect to the
%           coordinates of the cell. Since the hessian matrix is
%           symmetric only the upper half of feat is filled.
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
numFeatures = length(sigma)*(length(sigma) - 1)/2;
if isempty(raw)
    H = [];
    return;
end

nD = ndims(raw);
H = cell(nD,nD);
for i = 1:nD
    for j = i:nD
        order = zeros(ndims(raw),1);
        order(i) = order(i) + 1;
        order(j) = order(j) + 1;
        H{i,j} = Codat.Features.gaussianFilter(raw, sigma, truncate, order);
    end
end
end
