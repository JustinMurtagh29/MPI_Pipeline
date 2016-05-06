function grad = GaussGradient( raw, sigma, filterSize )
%GAUSSGRADIENTMAGNITUDE Calculate the gauss gradient of a 3D volume.
% INPUT sigma: Array containing standard deviation used for gaussian filter
%              in each dimension.
%       filterSize: Size of the resulting filter in each dimension and
%                   direction. The resulting filter has a size of
%                   2*filterSize + 1 and a total boundary of 2*filterSize.
% OUTPUT grad: cell array of length 3 containing the gradient for each
%              dimension
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

sigma = num2cell(sigma);
coords = arrayfun(@(siz)(-siz:siz).',filterSize,'UniformOutput',false);
coords = cellfun(@(coords,dim)permute(coords,[2:dim,1,(dim+1):3]),coords,num2cell(1:3),'UniformOutput',false);
gauss = cellfun(@(coords,sigma)exp(-(coords./sigma).^2./2),coords,sigma,'UniformOutput',false);
gauss = cellfun(@(gauss)gauss./sum(gauss),gauss,'UniformOutput',false);
gaussD = cellfun(@(coords,gauss,sigma)-gauss.*coords./(sigma.^2),coords,gauss,sigma,'UniformOutput',false);
grad = cellfun(@(gaussD)convn(raw,gaussD,'same'),gaussD,'UniformOutput',false);
for dim1 = 1:3
    for dim = setdiff(1:3,dim1)
        grad{dim1} = convn(grad{dim1},gauss{dim},'same');
    end
end


end

