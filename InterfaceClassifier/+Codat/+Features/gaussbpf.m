function [feat, border, numFeatures] = gaussbpf(raw,sigH,sigL)
%GAUSSBPF Perform gaussian band pass filtering in the fourier domain.
% Frequencies in the range [dh, dl] can pass the filter.
% INPUT raw: 3-dimensional input array.
%       sigH: [Nx3] or [Nx1] array of double specifying the standard
%           deviation for the gaussin high-pass filter in pixels for all
%           dimensions independently or jointly respectively for N
%           independent filterings of the raw data. If 0 is specified than
%           only low-pass filtering is performed (not optimized). N must be
%           the same for sigH and sigL.
%       sigL: [Nx3] or [Nx1] array of double specifying the standard
%           deviation for the gaussin low-pass filter in pixels for all
%           dimensions independently or jointly respectively. for N
%           independent filterings of the raw data. If Inf is specified
%           than only high-pass filtering is performed (not optimized). N 
%           must be the same for sigH and sigL.
% OUTPUT feat: [Nx1] cell array where N corresponds to the number of rows
%           in sigH and sigL. Each cell contains the filtered raw file with
%           the parameters specified in the corresponding rows of sigH and
%           sigL.
%        border: Integer array specifying the total border of the filter
%                in each dimension.
%        numFeatures: Resulting feature space size per voxel.
%
% Note The variance of the gaussians is chosen such that 2/3 of the signal
%      can pass at the specified standard deviation, i.e. the low-pass
%      filter has another factor 2 which is multiplied with the variance.
% 
% NOTE It is possible to calculate only the last two outputs by calling the
%      function without the raw input argument.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if size(sigH,1) ~= size(sigL,1)
    error('sigH and sigL must have the same number of rows.');
end

border = 0; %this is not correct since the filter has a boundary effect
numFeatures = size(sigH,1);
if isempty(raw)
    feat = [];
    return
end

if size(sigH,2) == 1
    sigH = repmat(sigH, [1, ndims(raw)]);
end
if size(sigL,2) == 1
    sigL = repmat(sigL, [1, ndims(raw)]);
end

%transform to fourier space
raw = fftn(raw);
raw = fftshift(raw);

%define grid for filter
[sx, sy, sz] = size(raw);
[x,y,z] = ndgrid(0:(sx - 1), 0:(sy - 1), 0:(sz - 1));
x = (x - floor(sx/2)).^2;
y = (y - floor(sy/2)).^2;
z = (z - floor(sz/2)).^2;

%prepare sig
sigH = sigH.^2;
sigL = 2.*sigL.^2;

feat = cell(numFeatures,1);
for i = 1:numFeatures
    %define filter
    lowPass = exp(-x./sigL(i,1) - y./sigL(i,2) - z./sigL(i,3));
    highPass = 1 - exp(-x./sigH(i,1) - y./sigH(i,2) - z./sigH(i,3));
    filter = lowPass.*highPass;
    clear lowPass highPass

    %filtering and backtransformation
    feat{i} = raw.*filter;
    feat{i} = ifftshift(feat{i});
    feat{i} = real(ifftn(feat{i}));
end


end
