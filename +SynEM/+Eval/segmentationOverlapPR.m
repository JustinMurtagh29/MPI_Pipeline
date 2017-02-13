function [ rp ] = segmentationOverlapPR( segGT, scores, thresholds, ...
    sizeTLower, sizeTUpper, overlapT, mergeTerm )
%SEGMENTATIONOVERLAPPR PR curve for segmentation overlap based performance.
% INPUT segGT: nd-array logical
%           Ground truth voxel labels.
%       scores: nd-array single
%           Voxelwise prediction scores/probabilities for ground truth
%           volume.
%       thresholds: [Nx1] float
%           Thresholds applied to scores.
%       sizeTLower: (Optional) double
%           see Interfaces.Eval.postProcessing
%           (Default: 0)
%       sizeTUpper: (Optional) double
%           see Interfaces.Eval.postProcessing
%           (Default: Inf)
%       overlapT: (Optional) double
%           Count overlaps only if they consist of overlapT voxels.
%           (Default: 0)
%       mergeTerm: (Optional) logical
%           Terminate if a predicted segment spans several ground truth
%           synapses.
%           (Default: true)
% OUTPUT rp: [Nx2] float
%           Recall-precision values for different thresholds.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('sizeTLower','var') || isempty(sizeTLower)
    sizeTLower = 0;
end
if ~exist('sizeTUpper','var') || isempty(sizeTUpper)
    sizeTUpper = Inf;
end
if ~exist('overlapT','var') || isempty(overlapT)
    overlapT = 0;
end
if ~exist('mergeTerm','var') || isempty(mergeTerm)
    mergeTerm = true;
end

segGTL = bwlabeln(segGT);
rp = cell(length(thresholds),1);
for i = 1:length(thresholds)
    currPred = SynEM.Eval.postProcessing(scores, thresholds(i), ...
        sizeTLower, sizeTUpper);
    [~,p, r] = SynEM.Eval.segmentationOverlapConfusionMat( ...
        segGTL, currPred, overlapT, mergeTerm);
    rp{i} = [r, p];
end
rp = cell2mat(rp);

end

