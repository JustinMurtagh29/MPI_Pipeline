function [binary_pred, CC, raw_crop] = postProcessing( prediction, prop_threshold, size_threshold, raw )
%POSTPROCESSING Post-processing on PSD probability/score map.
% Applies threshold on psd probability/score map and calculates connected
% components of output. Small connected components are discarded.
% INPUT prediction: PSD probability/score map.
%       prop_threshold: Threshold to apply on PSD probability/score map.
%       size_threshold: Minimal connected component size (in voxels).
%       raw: Optional raw input. Raw is cropped to size of prediction for
%            videos etc.
% OUTPUT binary_pred: Binary PSD prediction.
%        CC: Connected component information.
%        raw: Optional output if raw was specified.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

thresMap = prediction > prop_threshold;
CC = bwconncomp(thresMap,26);
CC = CC.PixelIdxList;
l = cellfun(@length,CC);
CC(l < size_threshold) = [];
CC = cell2struct(CC,'PixelIdxList',length(CC));

binary_pred = zeros(size(thresMap));
for i = 1:length(CC)
    binary_pred(CC(i).PixelIdxList) = 1;
    [x,y,z] = ind2sub(size(prediction),CC(i).PixelIdxList);
    CC(i).Centroid = [mean(x),mean(y),mean(z)];
end
binary_pred = logical(binary_pred);

if exist('raw','var') && nargout == 3
    border = (size(raw) - size(binary_pred))/2;
    raw_crop = raw(border(1) + 1:end - border(1),border(2) + 1:end - border(2),border(3) + 1:end - border(3));
elseif nargout == 3
    error('Input raw was not specified.\n');
end


end

