function [ Y ] = conv3fft( X, W, b, d, borderMode )
%CONV3FFT Convolution of video data with several channels by a bank of
% filters and adding an additional bias.
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
% see also conv3fft_2, conv3
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

sX = size(X);
sX(end+1:4) = 1;
w = flip(W,4); %compatibilty with convn version
sW = size(w);
Y = zeros([sX(1:3) - sW(1:3) + 1, sW(5)], 'like', X);
if prod(sW(1:3)) == 1  && numel(X)*sW(5)*4 < 2*1073741824
    Y = bsxfun(@plus,squeeze(sum(bsxfun(@times,X,w),permute(b,[4 3 2 1]))));
elseif prod(sW(1:3)) == 1
    for fm = 1:sW(5)
        Y(:,:,:,fm) = sum(bsxfun(@times,X,w(:,:,:,:,fm)),4) + b(fm);
    end
else % convolutional layer
    %determine size for fft (should be multiple of
    %sparsityFactor in all dimensions)
    
    filterSize = size(w);
    filterSize = filterSize(1:3);
    wSize = filterSize + (filterSize - 1).*(d - 1);
    
    switch borderMode
        case 'valid'
            ffSize = sX(1:3) + mod(-sX(1:3),d);
        case 'same'
            ffSize = sX(1:3) + ceil((wSize(1:3) - 1)./2) + mod(-(sX(1:3) + ceil((wSize(1:3) - 1)./2)),d);
        case 'full'
            ffSize = sX(1:3) + wSize(1:3) - 1 + mod(-(sX(1:3) + wSize(1:3) - 1),d);
    end
    
    %only do fft if filter size > 1
    fftDims = filterSize > 1;
    ffSize(~fftDims) = sX(~fftDims);
    fftDims = find(fftDims);
    for k = fftDims
        X = fft(X,ffSize(k),k);
        w = Codat.Lib.fftd(w,ffSize(k),k, d(k));
    end

    %fully vectorized product in fourier space
    %C = bsxfun(@times,activity,W);
    %C = permute(sum(C,4),[1, 2, 3, 5, 4]);

    %semi-vectorized with usually much less memory requirement and
    %almost same speed
    C = zeros([ffSize(1:3), sW(5)],'like',X);
    for i = 1:sW(5)
        C(:,:,:,i) = sum(bsxfun(@times,X,w(:,:,:,:,i)),4);
    end

    for k = fftDims
        C = ifft(C,ffSize(k),k);
    end

    if ~exist('borderMode','var') || isempty('borderMode')
        borderMode = 'valid';
    end

    %border behavior and bias
    switch borderMode
        case 'valid'
            Y = bsxfun(@plus,real(C(wSize(1):sX(1), ...
                wSize(2):sX(2),wSize(3):sX(3),:)), ...
                permute(b,[4 3 2 1]));
        case 'same'
            Y = bsxfun(@plus,real(C(ceil((wSize(1) - 1)/2) + (1:sX(1)), ...
                ceil((wSize(2) - 1)/2) + (1:sX(2)), ...
                ceil((wSize(3) - 1)/2) + (1:sX(3)),:)), ...
                permute(b,[4 3 2 1]));
        case 'full'
            Y = bsxfun(@plus,real(C), permute(b,[4 3 2 1]));
        otherwise
            error('Unknown border mode %s', borderMode);
    end

end


end
