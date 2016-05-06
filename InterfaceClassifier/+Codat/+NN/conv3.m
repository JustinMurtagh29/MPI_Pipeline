function [ Y ] = conv3( X, W, b, borderMode )
%CONV3 Convolution of video data with several channels by a bank of
% filters.
% INPUT X: 4d numerical array. The first three diensions correspond to the data
%           and the fourth to feature maps.
%       W: 5d numerical array. The first three dimensions correspond to the 3d
%           filter, the fourth dimension to the number of previous feature maps
%           and the fifth dimension to the number of output feature maps.
%       b: [Nx1] numerical array containing the bias for each output feature
%           map.
%       borderMode: (Optional) String specifying the convolution border
%           behaviour. Options are 'valid' (Default), 'same' and 'full'.
% OUTPUT Y: 4d numerical array. The first three dimensions correspond to the
%           data and the fourth to the output feature maps.
%
% see also conv3fft, conv3fft_2
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

sW = size(W);
sX = size(X);
switch borderMode
    case 'valid'
        Y = zeros([sX(1:3) - sW(1:3) + 1,sW(5)],'like',X);
        for fm = 1:sW(5)
            Y(:,:,:,fm) = convn(X,W(:,:,:,:,fm),borderMode) + b(fm);
        end
    case 'same'
        Y = zeros([sX(1:3),sW(5)],'like',X);
        for fm = 1:sW(5)
             tmp = convn(X,W(:,:,:,:,fm),borderMode) + b(fm);
             Y(:,:,:,fm) = tmp(:,:,:,ceil((sX(4) - 1)/2));
        end
    case 'full'
        Y = zeros([sX(1:3) + sW(1:3) - 1,sW(5)],'like',X);
        for fm = 1:sW(5)
             tmp = convn(X,W(:,:,:,:,fm),borderMode) + b(fm);
             Y(:,:,:,fm) = tmp(:,:,:,sX(4));
        end
end



end

