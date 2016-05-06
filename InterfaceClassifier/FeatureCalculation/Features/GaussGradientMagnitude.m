function imfeat = GaussGradientMagnitude( raw, sigma, filterSize )
%GAUSSGRADIENTMAGNITUDE Calculate the gauss gradient magnitude of a 3D volume.
% INPUT sigma: Array containing standard deviation used for gaussian filter
%              in each dimension.
%       filterSize: Size of the resulting filter in each direction. The
%                   resulting filter has a size of 2*filterSize + 1 and a
%                   total boundary of 2*filterSize.

grad = GaussGradient(raw,sigma,filterSize);
imfeat = sqrt(grad{1}.^2 + grad{2}.^2 + grad{3}.^2);

end

