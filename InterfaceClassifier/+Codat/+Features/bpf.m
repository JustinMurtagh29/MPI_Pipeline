function [ feat, border, numFeatures ] = bpf( raw, voxelSize, structureSize )
%BPF Band-pass filtering on 3d data using a gaussian band-pass filter in
%frequency space by converting the structure size to a frequency and
%applying the gaussianbpf function.
% INPUT raw: 3D input array.
%       voxelSize: [3x1] array of double specifying the physical size of a
%           voxel. Use [1, 1, 1] to specify the size in pixel.
%       structureSize: [Nx2] array of double specifying the lower and upper
%           in size for the target structure in physical units of voxel
%           size ranging from 0 to Inf for N independent filtering
%           operations.
% OUTPUT feat: [Nx1] cell array where N corresponds to the number of rows
%           in structureSize. Each cell contains an array of the same size
%           as raw with the filtered raw data.
%        border: Integer array specifying the total border of the filter
%                in each dimension.
%        numFeatures: Resulting feature space size per voxel.
%
% see also gaussbpf
% 
% NOTE It is possible to calculate only the last two outputs by calling the
%      function without the raw input argument.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%wavelength of structures (structure size is wavelength/2)
lowerWV = 2.*bsxfun(@rdivide,structureSize(:,1),voxelSize);
upperWV = 2.*bsxfun(@rdivide,structureSize(:,2),voxelSize);

%additional outputs
border = ceil(upperWV./2);
numFeatures = size(structureSize,2);
if isempty(raw)
    feat = [];
    return
end

%pixel frequency of structures
lowerF = 1./upperWV;
upperF = 1./lowerWV;

%frequency increment per pixel in fourier space
dF = 1./size(raw);

%region in fft space
sigH = bsxfun(@rdivide,lowerF,dF);
sigL = bsxfun(@rdivide,upperF,dF);

feat = Codat.Features.gaussbpf(raw, sigH, sigL);
end

