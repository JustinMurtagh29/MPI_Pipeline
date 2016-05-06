function S = structureTensor(raw, sigmaD, sigmaW, truncateD, truncateW)
%STRUCTURETENSOR Multidimensional structure tensor.
% The structure tensor is calculated for a gaussian derivative and a
% gaussian window function.
% INPUT raw: N-dimensional input array.
%       sigmaD: Integer scalar or vector specifying the standard deviation
%               of the gaussian derivative for all or seperately for each
%               dimension of the input image.
%       sigmaW: Integer scalar or vector specifying the standard deviation
%               of the gaussian window for all or seperately for each
%               dimension of the input image.
%       truncateD: (Optional) Double. Truncate filter at this many standard
%                 deviations.
%                 (Default: 3)
%       truncateW: (Optional) Double. Truncate filter at this many standard
%                 deviations.
%                 (Default: 3)
% OUTPUT S: Square cell array of size ndims(raw) in both dimension. Each
%           cell contains the second derivative with respect to the
%           coordinates of the cell. Since the struture tensor is
%           symmetric only the upper half of S is filled.
%
% NOTE This function uses convolutions in 'same' mode. The border for valid
%      convolutions for this filter is
%      (ceil(sigmaW*truncateW) - 1)./2 + (ceil(sigmaD*truncateD) - 1)./2.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%parse inputs
if ~isempty(raw) && length(sigmaD) == 1
    sigmaD = repmat(sigmaD,[1, ndims(raw)]);
elseif length(sigmaD) ~= ndims(raw)
    error('Specify the standard deviation for the derivative for all input dimensions.');
end
if ~exist('truncateD','var') || isempty(truncateD)
    truncateD = 3;
end
if ~isempty(raw) && length(sigmaW) == 1
    sigmaW = repmat(sigmaW,[1, ndims(raw)]);
elseif length(sigmaW) ~= ndims(raw)
    error('Specify the standard deviation for the window function for all input dimensions.');
end
if ~exist('truncateW','var') || isempty(truncateW)
    truncateW = 3;
end

nD = ndims(raw);
grad = cell(nD,1);
for dim = 1:ndims(raw)
    order = zeros(nD,1);
    order(dim) = 1;
    grad{dim} = Codat.Features.gaussianFilter(raw, sigmaD, truncateD, order);
end

S = cell(nD,nD);
for i = 1:nD
    for j = i:nD
        S{i,j} = Codat.Features.gaussianFilter(grad{i}.*grad{j}, sigmaW, truncateW);
    end
end
end
