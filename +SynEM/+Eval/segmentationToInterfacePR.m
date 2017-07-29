function [ rp ] = segmentationToInterfacePR( predMap, borders, ...
    groundTruth, thresholds, sizeT )
%SEGMENTATIONTOINTERFACEPR Calculate the PR curve of voxel-based prediction
%by determining synaptic interfaces using the voxel prediction map.
% INPUT predMap: 3d float
%       Voxel-wise PSD scores/probabilities.
%       borders: [Nx1] struct
%           Struct containing the borders for the predMap bounding box.
%       groundTruth: [Nx1] logical
%           Ground truth labels for borders.
%       thresholds: (Optional) [Nx1] float
%           Thresholds for which the rp points are calculated.
%           (Default: min(predMap(:)):0.01:max(predMap(:)))
%       sizeT: (Optional) double
%           Minimal connected component size of predicted components.
%           (Default: 50)
% OUTPUT rp: [Nx2] float
%           Recall-precision values.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('thresholds','var') || isempty(thresholds)
    thresholds = min(predMap(:)):0.01:max(predMap(:));
end
if ~exist('sizeT','var') || isempty(sizeT)
    sizeT = 50;
end

rp = zeros(length(thresholds),1);
warning('off','getSynapticBorders:w1');
tic
for i = 1:length(thresholds)
    tmp = SynEM.Eval.postProcessing(predMap, thresholds(i), sizeT);
    % y = SynEM.Eval.getSynapticBorders(borders, tmp, 0.05);
    y = SynEM.Eval.getSynapticBordersMaxOverlap(borders, tmp);
    [~,rp(i,2), rp(i,1)] = SynEM.Util.confusionMatrix(groundTruth, y);
    Util.progressBar(i, length(thresholds));
end

end
