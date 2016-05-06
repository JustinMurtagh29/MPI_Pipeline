function [ raw, target, targetWeights ] = loadDataMT( path )
%LOADDATAMT Load data for multiple targets and tanh outputs.
% INPUT path: Path to training data file.
% REMINDER targets labels (i.e. file in path):
%                 0: intracellular
%                 1: membrane
%                 2: vesicles
%                 3: mitos
% REMINDER target labels (i.e. output target)
%                 -1: membrane
%                 1: intracellular
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

load(path);
raw = (single(raw) - 122)./22;
tmp_target = zeros([size(targets),3],'single');
for i = 1:3
    labels = targets == i;
    tmp_target(:,:,:,i) = labels - ~labels;
end
targetWeights = ~repmat(target == 0,[1 1 1 3]);
target = tmp_target;

end

