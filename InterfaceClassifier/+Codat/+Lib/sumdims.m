function [ X ] = sumdims( X, dims )
%SUMDIMS Like Matlab sum for several dimensions.
% INPUT X: Numerical array of arbitrary shape.
%       dims: [Nx1] array of integer specifying the dimensions for which
%           the sum is calculated.
% OUTPUT X: The input array summed over dims.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

for dim = dims
    X = sum(X,dim);
end
end

