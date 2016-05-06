function [ raw, target, targetWeights ] = loadDataMTMasked( path )
%LOADDATAMTMASKED Load data for multiple targets and tanh outputs. If the
% datafile does not contain vesicle and mito labels they are masked.
% INPUT path: Path to training data file.
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

load(path);
raw = (single(raw) - 122)./22;
if exist('targets','var')
    tmp_target = zeros([size(targets),3],'single');
    for i = 1:3
        labels = targets == i;
        tmp_target(:,:,:,i) = labels - ~labels;
    end
    targetWeights = repmat(target == 0,[1 1 1 3]);
    target = tmp_target;
else
    targetWeights = false([size(target),3]);
    targetWeights(:,:,:,1) = target == 0;
    targetWeights(:,:,:,2:3) = true;
    target = -single(target);
    target(:,:,:,2:3) = zeros([size(target),2]);
end
targetWeights = ~targetWeights;

end

