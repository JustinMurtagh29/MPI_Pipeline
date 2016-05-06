function [feat, border, numFeatures] = eigenvaluesOfStructureTensor( raw, sigmaD, sigmaW, truncateD, truncateW)
%EIGENVALUESOFSTRUCTURETENSOR Calculate the eigenvalues of a 3x3 structure
% tensor.
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
% OUTPUT feat: Cell array of length 3 containing the eigenvalues in
%              increasing order.
%        border: Integer array specifying the total border of the filter
%                in each dimension.
%        numFeatures: Resulting feature space size per voxel.
%
% NOTE It is possible to calculate only the last two outputs by calling the
%      function without the raw input argument.
% NOTE This function uses convolutions in 'same' mode. The border for valid
%      convolutions for this filter is
%      (ceil(sigmaW*truncateW) - 1)./2 + (ceil(sigmaD*truncateD) - 1)./2.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%parse inputs
if ~isempty(raw) && length(sigmaD) == 1
    sigmaD = repmat(sigmaD,[1, ndims(raw)]);
elseif ~isempty(raw) && length(sigmaD) ~= ndims(raw)
    error('Specify the standard deviation for the derivative for all input dimensions.');
end
if ~exist('truncateD','var') || isempty(truncateD)
    truncateD = 3;
end
if ~isempty(raw) && length(sigmaW) == 1
    sigmaW = repmat(sigmaW,[1, ndims(raw)]);
elseif ~isempty(raw) && length(sigmaW) ~= ndims(raw)
    error('Specify the standard deviation for the window function for all input dimensions.');
end
if ~exist('truncateW','var') || isempty(truncateW)
    truncateW = 3;
end
border = 2.*(ceil(sigmaW*truncateW)  + ceil(sigmaD*truncateD));
numFeatures = 3;
if isempty(raw)
    feat = [];
    return;
end

S = Codat.Features.structureTensor(raw, sigmaD, sigmaW, truncateD, truncateW);

[nx, ny, nz] = size(raw);
eigenvalues = eig3S(S{1,1}(:),S{1,2}(:),S{1,3}(:),S{2,2}(:),S{2,3}(:),S{3,3}(:));
feat = cell(3,1);
feat{3} = reshape(eigenvalues(:,1),nx,ny,nz);
feat{2} = reshape(eigenvalues(:,2),nx,ny,nz);
feat{1} = reshape(eigenvalues(:,3),nx,ny,nz);
end
