function [raw, target, targetWeights] = loadDataRep(path)
%LOADDATAREP Data for tanh output layer with input data replication.
% INPUT path: Path to training data file.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

m = load(path);
raw = repmat((single(m.raw) - 122)./22,[1 1 1 8]);
target = single(m.target);
targetWeights = ones(size(target),'single');
targetWeights(target == -1) = 3;
targetWeights(target == 0) = 0;
end
