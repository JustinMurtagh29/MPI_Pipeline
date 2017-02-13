function [seg, stats, raw_crop] = postProcessing( prediction, probT, ...
    sizeTLower, sizeTUpper, raw )
%POSTPROCESSING Post-processing on PSD probability/score map.
% Applies threshold on psd probability/score map and calculates connected
% components of output. Small connected components are discarded.
% INPUT prediction: 3d float
%           Prediction probability/score map.
%       probT: double
%           Threshold to apply on PSD probability/score map.
%       sizeTLower: int
%           Minimal connected component size (in voxels).
%       sizeTUpper: (Optional) int
%           Maximal connected  component size (in voxels).
%           (Default: Inf).
%       raw: (Optional) 3d numerical 
%           Raw is cropped to size of prediction for videos etc.
% OUTPUT seg: 3d int
%           Predicted segmentation.
%        stats: struct
%           Struct containing the connected component information.
%           Fields are 'Area','Centroid','PixelIdxList'.
%        raw: 3d numerical
%           Optional output if raw was specified.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('sizeTUpper','var') || isempty(sizeTUpper)
    sizeTUpper = Inf;
end

bw = prediction >= probT;
bw = imclose(bw,ones(3, 3, 3));
stats = regionprops(bw,'Area','Centroid','PixelIdxList');
stats([stats.Area] < sizeTLower) = [];
stats([stats.Area] > sizeTUpper) = [];
seg = false(size(bw));
seg(cell2mat({stats.PixelIdxList}')) = true;

if exist('raw','var') && nargout == 3
    border = (size(raw) - size(seg))/2;
    raw_crop = raw(border(1) + 1:end - border(1), ...
        border(2) + 1:end - border(2),border(3) + 1:end - border(3));
elseif nargout == 3
    warning('Third output is empty since raw input was not specified.\n');
    raw_crop = [];
end

end

