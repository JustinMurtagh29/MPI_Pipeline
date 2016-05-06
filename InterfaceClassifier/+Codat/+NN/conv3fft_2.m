function [ Y ] = conv3fft_2( X, W, b, d, borderMode )
%CONV3FFT_2 As conv3fft but iterated over feature maps for less memory
%requirement.
% INPUT X: 4d numerical array. The first three diensions correspond to the data
%           and the fourth to feature maps.
%       W: 5d numerical array. The first three dimensions correspond to the 3d
%           filter, the fourth dimension to the number of previous feature maps
%           and the fifth dimension to the number of output feature maps.
%       b: [Nx1] numerical array containing the bias for each output feature
%           map.
%       d: [3x1] array of integer containing the d-regularity factor for each
%           dimension.
%       borderMode: (Optional) String specifying the convolution border
%           behaviour. Options are 'valid' (Default), 'same' and 'full'.
% OUTPUT Y: 4d numerical array. The first three dimensions correspond to the
%           data and the fourth to the output feature maps.
%
% see also conv3fft, conv3
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%convolution via fft iterated over fms (medium memory
%requirement)
W = W(:,:,:,end:-1:1,:); %compatibility with conv3 version
sX = size(X);
sX(end+1:4) = 1;
sW = size(W);
filterSize = sW(1:3);
wSize = filterSize + (filterSize - 1).*(d - 1);

switch borderMode
    case 'valid'
        ffSize = sX(1:3) + mod(-sX(1:3),d);
        Y = zeros([sX(1:3) - wSize(1:3) + 1,sW(5)],'like',X);
    case 'same'
        ffSize = sX(1:3) + ceil((wSize(1:3) - 1)./2) + mod(-(sX(1:3) + ceil((wSize(1:3) - 1)./2)),d);
        Y = zeros([sX(1:3),sW(5)],'like',X);
    case 'full'
        ffSize = sX(1:3) + wSize(1:3) - 1 + mod(-(sX(1:3) + wSize(1:3) - 1),d);
        Y = zeros([sX(1:3) + wSize(1:3) - 1,sW(5)],'like',X);
end

%only do fft if filter size > 1
fftDims = filterSize > 1;
ffSize(~fftDims) = sX(~fftDims);
fftDims = find(fftDims);

%fft activity
for k = fftDims
    X = fft(X,ffSize(k),k);
end

%fft weights and ifft activity*weights for current feature map

for fm = 1:sW(5)
    w = W(:,:,:,:,fm);
    for k = fftDims
        w = Codat.Lib.fftd(w,ffSize(k),k, d(k));
    end
    C = sum(bsxfun(@times,X,w),4); %bsxfun necessary if any dimension in w is 1
    for k = fftDims
        C = ifft(C,ffSize(k),k);
    end
    
    %border behavior and bias
    switch borderMode
        case 'valid'
            Y(:,:,:,fm) = real(C(wSize(1):sX(1), ...
                wSize(2):sX(2),wSize(3):sX(3))) + b(fm);
        case 'same'
            Y(:,:,:,fm) = real(C(ceil((wSize(1) - 1)/2) + (1:sX(1)), ...
                ceil((wSize(2) - 1)/2) + (1:sX(2)), ...
                ceil((wSize(3) - 1)/2) + (1:sX(3)))) + b(fm);
        case 'full'
            Y(:,:,:,fm) = real(C) + b(fm);
        otherwise
            error('Unknown border mode %s', borderMode);
    end
end


end

