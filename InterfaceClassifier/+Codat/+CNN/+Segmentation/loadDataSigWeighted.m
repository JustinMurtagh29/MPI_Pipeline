function [ raw, target, targetWeights ] = loadDataSigWeighted( path, weight )
%LOADDATASIG Data for sigmoid output layer.
% INPUT path: Path to training data file.
%       weight: (Optional) Double weight factor for membran class.
%               (Default: Emprical weight to balance classes)
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

m = load(path);

if ~exist('weight','var') || isempty(weight)
    weight = sum(m.target(:) == 1)/sum(m.target(:) == -1);
end

raw = (single(m.raw) - 122)./22;
tmp = zeros(size(m.target),'single');
%m.target = imdilate(m.target,ones(2,2));
tmp(m.target == 1) = 0;
tmp(m.target == -1) = 1; %membrane
target = tmp;
targetWeights = ones(size(target),'single');
targetWeights(target == 1) = weight;
targetWeights(m.target == 0) = 0;


end

