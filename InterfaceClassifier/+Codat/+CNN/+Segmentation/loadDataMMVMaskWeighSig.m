function [ raw, target, targetWeights ] = loadDataMMVMaskWeighSig( path, weight )
%LOADDATAMMVMASKWEIGHSIG Load data for multiple targets and tanh outputs.
% If the datafile does not contain vesicle and mito labels they are masked.
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
if isfield(m,'targets')
    tmp_target = zeros([size(m.targets),3],'single');
    for i = 1:3
        tmp_target(:,:,:,i) = m.targets == i;
    end
    targetWeights = ones([size(m.targets),3],'single');
    targetWeights(:,:,:,1) = weight.*single(m.target == -1) + single(m.target == 1) + 0.*single(m.target == 0);
    target = tmp_target;
else
    targetWeights = zeros([size(m.target),3]);
    targetWeights(:,:,:,1) = weight.*single(m.target == -1) + single(m.target == 1);
    target = zeros([size(m.target),3],'single');
    target(m.target == -1) = 1;
end


end

