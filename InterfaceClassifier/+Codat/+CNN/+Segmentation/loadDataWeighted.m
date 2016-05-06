function [ raw, target, targetWeights ] = loadDataWeighted( path, weight, sT )
%LOADDATAWEIGHTED Data for tanh output layer.
% INPUT path: Path to training data file.
%       weight: (Optional) Double weight factor for membran class.
%               (Default: Emprical weight to balance classes)
%       sT: (Optional) Logical specifying whether to train on small target.
%           (Default: false)
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

m = load(path);

if ~exist('weight','var') || isempty(weight)
    weight = sum(m.target(:) == 1)/sum(m.target(:) == -1);
end

raw = (single(m.raw) - 122)./22;
target = single(m.target);
if exist('sT','var') && sT
    target = imdilate(target,ones(2,2));
end
targetWeights = ones(size(target),'single');
targetWeights(target == -1) = weight;
targetWeights(target == 0) = 0;

end
