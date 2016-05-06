function [T,precision, recall] = centroidBasedPerformance( CC, targetCentroids, raw, pred, distanceThreshold )
%CENTROIDBASEDPERFORMANCE Estimate performance of synapse detection based
% on distance of centroids of connected components of predictions and
% centroids of synaptic interfaces.
% INPUT CC: Connected components returned from postProcessing of pred.
%       targetCentroids: Centroids of targets in raw.
%       raw: Input raw data.
%       prediction: PSD probability/score map.
%       distanceThreshold: Maximal distance between predicted and ground
%           truth synapse in nm.
% OUTPUT T: confusion matrix (fn Nan because it can not be calculated).
%        precision: Precision calculated from T.
%        recall: Recall calculated from T.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%calculate CC centroids in coordinates of raw
syn_pred = reshape([CC(:).Centroid],3,length(CC))';
translation = (size(raw) - size(pred))/2;
syn_pred = bsxfun(@plus,syn_pred,translation);

%convert coordinates to nm
syn_pred = bsxfun(@times,syn_pred,[11.24 11.24 28]);
targetCentroids = bsxfun(@times,targetCentroids,[11.24 11.24 28]);

%calculate precision and recall based on max distance threshold
D = pdist2(syn_pred,targetCentroids);
pred = D < distanceThreshold;
T(1,1) = sum(any(pred,2),1); %tp
T(1,2) = sum(all(~pred,2),1); %fp
T(2,1) = sum(all(~pred,1),2); %fn
T(2,2) = NaN; %tn
precision = T(1,1)/(T(1,1) + T(1,2));
recall = T(1,1)/(T(1,1) + T(2,1));
end

