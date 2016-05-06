function [ raw, target, targetWeights ] = loadDataForMembranePredictionTest( path )
%LOADDATAFORMEMBRANEPREDICTIONTEST Test membrane prediction for MT cnn.
% INPUT path: Path to training data file.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

load(path);
raw = (single(raw) - 122)./22;
target = -single(target);
target = cat(4,target,zeros([size(target), 2]));
targetWeights = ~(target == 0);

end