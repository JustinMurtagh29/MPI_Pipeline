function [rp, thresholds] = pointDistanceDetectionRP( pred_coords, ...
    scores, gt_coords, maxDist, voxelSize, thresholds, bbox )
%POINTDISTANCEDETECTIONRP Synapse detection performance based on maximal
%distance of predicted to ground truth synapse coms.
% TPs are defined by a gt coordinate that has a pred coordinates within
% maxDist.
% FPs are defined by pred coordinates that do not have a gt coordinate
% within maxDist.
% FNs are defined by gt coordinates that do not have a gt coordinates
% within maxDist.
% INPUT pred_coords: [Nx3] int
%           The coordinates of predicted synapses.
%       scores: [Nx1] float
%           The scores for each pred_coord.
%       gt_coords: [Nx3] int
%           Coordinates of ground truth synapses.
%       maxDist: double
%           Maximal distance in the units of voxel size.
%       voxelSize: (Optional) [1x3] double
%           Size of a voxel in each dimension.
%           (Default: [1, 1, 1])
%       thresholds: (Optional) [Nx1] float
%           Thre thresholds to evaluate.
%           (Default: min(unique(scores)):0.01:(max(unique(scores)) + 0.01)
%       bbox: (Optional) [3x2] int
%           Restrict to synapses within the specified bounding box. GTs
%           outside of the bounding box do not cause FNs or FPs, i.e. they
%           need not be found and detections within maxDist are not counted
%           as FPs.
%           (Default: no restriction of gt nodes.)
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

assert(size(pred_coords, 1) == length(scores), ...
    'Specify a score for each prediction.');

% get gt nodes out of bbox (oob)
if exist('bbox', 'var') && ~isempty(bbox)
    oob_gt = ~Util.isInBBox(gt_coords, bbox);
    oob_pred = ~Util.isInBBox(pred_coords, bbox);
else
    oob_gt = false(size(gt_coords, 1), 1);
    oob_pred = false(size(pred_coords, 1), 1);
end

% scale coordinates by voxel size
if exist('voxelSize', 'var') && ~isempty(voxelSize)
    voxelSize = voxelSize(:)';
    pred_coords = bsxfun(@times, double(pred_coords), voxelSize);
    gt_coords = bsxfun(@times, double(gt_coords), voxelSize);
end

if ~exist('thresholds', 'var') || isempty(thresholds)
    thresholds = unique(scores);
    thresholds = min(thresholds):0.01:(max(thresholds) + 0.01);
end
    
d = pdist2(pred_coords, gt_coords) <= maxDist;

tps = zeros(length(thresholds), 1);
fns = zeros(length(thresholds), 1);
fps = zeros(length(thresholds), 1);
for i = 1:length(thresholds)
    cur_preds = scores >= thresholds(i);
    curD = d(cur_preds, :);
    tps(i) = sum(any(curD(:,~oob_gt), 1));
    fns(i) = sum(~any(curD(:,~oob_gt), 1));
    fps(i) = sum(~any(curD(oob_pred(cur_preds),:), 2));
end

rp = [tps./(tps + fns), tps./(tps + fps)];

end

