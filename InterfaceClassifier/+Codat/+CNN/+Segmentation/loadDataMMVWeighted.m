function [ raw, target, targetWeights ] = loadDataMMVWeighted( path, weight )
%LOADDATAMMVWEIGHTED Load MMV data for tanh output.
% Membrane labels are weighted with a factor.
% INPUT path: Path to training data file.
%       weight: (Optional) Double weight factor for membran class.
%               (Default: Emprical weight to balance classes)
% REMINDER targets Labels:
%                 0: intracellular
%                 1: membrane
%                 2: vesicles
%                 3: mitos
% REMINDER target labels:
%                 1: intracellular
%                 0: mask
%                 -1: membrane
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

m = load(path);
if ~exist('weight','var') || isempty(weight)
    weight = sum(m.target(:) == 1)/sum(m.target(:) == -1);
end

raw = (single(m.raw) - 122)./22;
tmp_target = zeros([size(m.targets),3],'single');
for i = 1:3
    tmp_target(:,:,:,i) = (m.targets == i) - ~(m.targets == i);
end
targetWeights = ones([size(m.targets),3],'single');
targetWeights(:,:,:,1) = weight.*single(m.target == -1) + single(m.target == 1) + 0.*single(m.target == 0);
target = tmp_target;

end

