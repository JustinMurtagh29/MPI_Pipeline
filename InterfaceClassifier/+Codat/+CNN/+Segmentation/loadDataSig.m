function [ raw, target, targetWeights ] = loadDataSig( path )
%LOADDATASIG Data for sigmoid output layer.
% INPUT path: Path to training data file.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

m= load(path);
raw = (single(m.raw) - 122)./22;
tmp = zeros(size(m.target),'single');
tmp(m.target == -1) = 1; %membrane
tmp(m.target == 1) = 0;
targetWeights = ~(m.target == 0);
target = tmp;


end

