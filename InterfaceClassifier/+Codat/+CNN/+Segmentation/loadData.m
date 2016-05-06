function [ raw, target, targetWeights ] = loadData( path )
%LOADDATA Data for tanh output layer.
% INPUT path: Path to training data file.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

load(path);
raw = (single(raw) - 122)./22;
target = single(target);
targetWeights = ~(target == 0);

end

