function [ raw, target, targetWeights ] = loadDataWeighted( path )
%LOADDATAWEIGHTED Data for tanh output layer where area around psd is
% masked and weighted to balance classes.
% INPUT path: Path to training data file.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

load(path);
raw = (single(raw) - 122)./22;
target = single(target);

targetMask = false(size(target));
for i = 1:size(target,3)
    targetMask(:,:,i) = imdilate(target(:,:,i),ball(6,[1 1 2.5]));
end
targetMask(target > 0) = false;
synFactor = numel(target)/sum(target(:));
targetWeights = ones(size(target),'single');
targetWeights(target == 1) = synFactor;
targetWeights(targetMask) = 0;
normalization = sum(targetMask(:))./sum(targetWeights(:));
targetWeights = targetWeights.*normalization;
target(target == 0) = -1;

end

function h = ball( radius, factor )
%BALL
[x,y,z] = meshgrid(-ceil(radius/factor(1)):ceil(radius/factor(1)), ...
                   -ceil(radius/factor(2)):ceil(radius/factor(2)), ...
                   -ceil(radius/factor(3)):ceil(radius/factor(3)));
h = (factor(1).*x).^2 + (factor(2).*y).^2 + (factor(3).*z).^2 <= radius^2;
end
