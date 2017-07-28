function [rp, thresholds] = segmentationPointDistanceDetectionRP( ...
    pred, gt_coords, maxDist, voxelSize, thresholds, bbox  )
%SEGMENTATIONPOINTDISTANCEDETECTIONRP Synapse detection performance based
% on maximal distance of predicted to ground truth synapse coms. Coms are
% extracted from the threshold input prediction.
% TPs are defined by a gt coordinate that has a pred coordinates within
% maxDist.
% FPs are defined by pred coordinates that do not have a gt coordinate
% within maxDist.
% FNs are defined by gt coordinates that do not have a gt coordinates
% within maxDist.
% INPUT pred: 3d float
%           The voxel-wise prediction map.
%       gt_coords: [Nx3] int
%           Coordinates of ground truth synapses.
%       maxDist: double
%           Maximal distance in the units of voxel size.
%       voxelSize: (Optional) [1x3] double
%           Size of a voxel in each dimension.
%           (Default: [1, 1, 1])
%       thresholds: (Optional) [Nx1] float
%           Thre thresholds to evaluate.
%           (Default: (min(pred(:)) + 0.05):0.05:max(pred(:)))
%       bbox: (Optional) [3x2] int
%           Restrict to synapses within the specified bounding box. GTs
%           outside of the bounding box do not cause FNs or FPs, i.e. they
%           need not be found and detections within maxDist are not counted
%           as FPs.
%           (Default: no restriction of gt nodes.)
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

% get gt nodes out of bbox (oob)
if exist('bbox', 'var') && ~isempty(bbox)
    oob = ~Util.isInBBox(gt_coords, bbox);
else
    oob = false(size(gt_coords), 1);
end

% scale coordinates by voxel size
if exist('voxelSize', 'var') && ~isempty(voxelSize)
    voxelSize = voxelSize(:)';
    gt_coords = bsxfun(@times, double(gt_coords), voxelSize);
else
    voxelSize = [];
end

if ~exist('thresholds', 'var') || isempty(thresholds)
    thresholds = (min(pred(:))+ 0.05):0.05:max(pred(:));
end

tps = zeros(length(thresholds), 1);
fns = zeros(length(thresholds), 1);
fps = zeros(length(thresholds), 1);
for i = 1:length(thresholds)
    pred_ps = SynEM.Eval.postProcessing(pred, thresholds(i), 50);
    stats = regionprops(pred_ps, 'Centroid');
    pred_coords = round(reshape([stats.Centroid], 3, [])');
    pred_coords = pred_coords(:,[2 1 3]);
    if ~isempty(voxelSize)
        pred_coords = bsxfun(@times, double(pred_coords), voxelSize);
    end
    curD = pdist2(pred_coords, gt_coords) <= maxDist;
    tps(i) = sum(any(curD(:,~oob), 1));
    fns(i) = sum(~any(curD(:,~oob), 1));
    fps(i) = sum(~any(curD, 2));
end

rp = [tps./(tps + fns), tps./(tps + fps)];

end

