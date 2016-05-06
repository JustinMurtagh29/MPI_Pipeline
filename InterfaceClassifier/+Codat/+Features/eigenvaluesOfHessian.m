function [feat, border, numFeatures] = eigenvaluesOfHessian(raw, sigma, truncate)
%EIGENVALUESOFHESSIAN Calculate the eigenvalues of a 3x3 hessian matrix.
% INPUT raw: 3-dimensional input array.
%       sigma: Integer scalar or vector specifying the standard deviation of
%              the gaussian kernel for all or seperately for each dimension
%              of the input image.
%       truncate: (Optional) Double. Truncate filter at this many standard
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
numFeatures = 3;
if isempty(raw)
    feat = [];
    return;
end

H = Codat.Features.hessianMatrix(raw, sigma, truncate);

[nx, ny, nz] = size(raw);
eigenvalues = eig3S(H{1,1}(:),H{1,2}(:),H{1,3}(:),H{2,2}(:),H{2,3}(:),H{3,3}(:));
feat = cell(3,1);
feat{3} = reshape(eigenvalues(:,1),nx,ny,nz);
feat{2} = reshape(eigenvalues(:,2),nx,ny,nz);
feat{1} = reshape(eigenvalues(:,3),nx,ny,nz);
end
