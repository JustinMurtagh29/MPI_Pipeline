function imfeat = DifferenceOfGaussians( raw, sigma, filterSize, k )
%DIFFERENCEOFGAUSSIANS Calculate DoG of a 3D volume (image stack).
% INPUT sigma: Array containing standard deviation used for gaussian filter
%              in each dimension.
%       filterSize: Size of the resulting filter in each dimension and
%                   direction. The resulting filter has a size of
%                   2*filterSize + 1 and a total boundary of 2*filterSize.
%       k: Factor 1/k in standard deviation of second gaussian (Must be > 0).

sigma = num2cell(sigma);
coords = arrayfun(@(siz)(-siz:siz).',filterSize,'UniformOutput',false);
coords = cellfun(@(coords,dim)permute(coords,[2:dim,1,(dim+1):3]),coords,num2cell(1:3),'UniformOutput',false);
gauss1 = cellfun(@(coords,sigma)exp(-(coords./sigma).^2./2),coords,sigma,'UniformOutput',false);
gauss1 = cellfun(@(gauss)gauss./sum(gauss),gauss1,'UniformOutput',false);
gauss2 = cellfun(@(coords,sigma)exp(-(coords./sigma).^2./2./k^2),coords,sigma,'UniformOutput',false);
gauss2 = cellfun(@(gauss)gauss./sum(gauss),gauss2,'UniformOutput',false);
I1 = raw;
for dim = 1:3
    I1 = convn(I1,gauss1{dim},'same');
end
I2 = raw;
for dim = 1:3
    I2 = convn(I2,gauss2{dim},'same');
end
imfeat = I1 - I2;


end

