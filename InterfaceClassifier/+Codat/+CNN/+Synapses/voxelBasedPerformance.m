function [T,recall] = voxelBasedPerformance( pred, targetCentroids, distance, numVoxels )
%VOXELBASEDPERFORMANCE Estimate recall of synapse detection based on
%psd voxels close to centroids of ground truth synapses. This should be
%used for performance evaluation with very small thresholds.
% INPUT pred: Binary psd prediction map.
%       targetCentroids: Centroids of ground truth synapses as n x 3
%           matrix.
%       distance: Distance around target centroids to look for psd
%           prediction in each dimension.
%       numVoxels: Number of voxels within distance around target centroids
%           to accept target as found.
% OUTPUT recall:

targetCentroids = round(targetCentroids);
wasFound = false(length(targetCentroids),1);
for i = 1:length(targetCentroids)
    cen = targetCentroids(i,:);
    probs = pred(max(cen(1) - distance(1),1): min(cen(1) + distance(1),size(pred,1)), ...
                 max(cen(2) - distance(2),1): min(cen(2) + distance(2),size(pred,2)), ...
                 max(cen(3) - distance(3),1): min(cen(3) + distance(3),size(pred,3)));
    wasFound(i) = sum(probs(:)) >= numVoxels;
end
T(1,1) = sum(wasFound);
T(2,1) = sum(~wasFound);
recall = sum(wasFound)/length(wasFound);

end

