function [ T, p, r ] = segmentationOverlapConfusionMat( segGT, pred, ...
    overlapT, mergeTerm, nhood )
%SEGMENTATIONOVERLAPCONFUSIONMAT Calculate the performance of a classifier
%by determining the overlap of ground truth segmentation and predicted
%segmentation.
% INPUT segGT: nd array logical or int
%           Logical array corresponding to the ground truth segmentation
%           or ground truth as label matrix.
%       pred: nd-array logical or int
%           Logical array containing the predicted segmentation of same
%           size .as ground truth or a label matrix.
%       overlapT: (Optional) double
%           Count overlaps only if they consist of overlapT voxels.
%           (Default: 0)
%       mergeTerm: (Optional) logical
%           Terminate if a predicted segment spans several ground truth
%           synapses.
%           (Default: true)
%       nhood: (Optional) [1x3] int
%           Neighborhood for image dilation of segGT.
%           (Default: no dilation).
% OUTPUT T: Confusion matrix.
%        p: Precision.
%        r: Recall.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('overlapT','var') || isempty(overlapT)
    overlapT = 0;
end

if ~exist('mergeTerm','var') || isempty(mergeTerm)
    mergeTerm = true;
end

if exist('nhood','var') && ~isempty(nhood)
    if iscolumn(nhood)
        nhood = nhood';
    end
    segGT = imdilate(segGT, ones(nhood));
end

if islogical(segGT)
    [segGTL, numGT] = bwlabeln(segGT);
else
    segGTL = segGT;
    segGT = segGTL > 0; %not really necessary
    numGT = max(segGTL(:));
end

if islogical(pred)
    detectedSynapticPixels = segGT & pred;
    [predL,numPred] = bwlabeln(pred);
else
    detectedSynapticPixels = segGT & (pred > 0);
    predL = pred;
    numPred = max(predL(:));
end

if mergeTerm
    %check if a predicted cluster has overlap with several GT clusters
    stats = regionprops(predL,'PixelIdxList');
    stats = {stats.PixelIdxList};
    numGTClusters = cellfun(@(x)length(unique(setdiff(segGTL(x),0))), ...
                            stats);
   if any(numGTClusters > 1) %terminate
       T = nan(2,2);
       p = nan;
       r = nan;
       return
   end
end


%only count a detection if at least overlapT pixels overlap
detPixels = predL(detectedSynapticPixels);
detClusters = unique(detPixels);
n = histc(detPixels,detClusters);
detClusters = detClusters(n >= overlapT);

%same for ground truth clusters
detPixels = segGTL(detectedSynapticPixels);
detClustersGT = unique(detPixels);
n = histc(detPixels,detClustersGT);
detClustersGT = detClustersGT(n >= overlapT);

%calculate confusion matrix entries
tpClusters = detClustersGT;
fpClusters = setdiff(1:numPred,detClusters);
% tpClustersGT = setdiff(detClustersGT,0);
% tpClustersGT = setdiff(unique(segGTL(detectedSynapticPixels)),0);
fnClusters = setdiff(1:numGT,tpClusters);

%outputs
T = [length(tpClusters), length(fpClusters); length(fnClusters), NaN];
p = T(1,1)/(T(1,1) + T(1,2));
r = T(1,1)/(T(1,1) + T(2,1));

end
