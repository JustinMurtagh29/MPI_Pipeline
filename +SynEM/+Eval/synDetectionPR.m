function rp = synDetectionPR( segGT, synScores, borders )
%SYNDETECTIONPR Synapse detection performance on ground truth segmentation.
% INPUT segGT: Logical array corresponding to the ground truth segmentation.
%       synScores: [Nx1] float array containing the synapse scores for each
%           interface contained in the segGT bounding box.
%       borders: [Nx1] struct containing the field PixelIdxList with linear
%           indices of border voxels w.r.t. to segGT.
% OUTPUT rp: [Nx2] float containing the recall and precision values for the
%           test set.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

thresholds = min(synScores):0.05:max(synScores);
rp = cell(length(thresholds),1);
segGTL = bwlabeln(segGT);
parfor i = 1:length(thresholds)
    synPred = synScores >= thresholds(i);
    predMat = false(size(segGT));
    predMat(cell2mat({borders(synPred).PixelIdxList})) = true;
    [~,p,r] = SynEM.Eval.segmentationOverlapPerformance(segGTL, predMat);
    rp{i} = [r, p];
end
rp = cell2mat(rp);

end

