function h = plotWeights( cnet, layer, featureMaps)
%PLOTWEIGHTS Basic plotting of cnet weights.
% INPUT layer: Specify layer.
%       featureMaps: Specify which feature map to plot. (Default: all feature
%                 maps).
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if layer < 2 || layer > cnet.layer
    h = [];
    return;
end

if ~exist('featureMaps','var') || isempty(featureMaps)
    featureMaps = 1:cnet.featureMaps(layer);
end

weightImage = [];
for fm = featureMaps
    column = [];
    for prevFm = 1:cnet.featureMaps(layer - 1)
        row  = zeros(length(cnet.W{layer}(:,1,1,1,fm)),2);
        for frame = 1:size(cnet.W{layer}(:,:,:,prevFm,fm),3)
            row = cat(2,row,cnet.W{layer}(:,:,frame,prevFm,fm),zeros(length(cnet.W{layer}(:,1,frame,prevFm,fm)),2));
        end
        column = cat(1,column,zeros(2,size(row,2)),row);
    end
    weightImage = cat(1,weightImage,column,zeros(2,size(row,2)));
end
h = figure;
imagesc(weightImage);
title(sprintf('Weight matrices - layer %d',layer));
set(gca,'XTickLabel','Filter z-dimension');
set(gca,'XTick',size(weightImage,2)/2);
set(gca,'YTickLabel','Selected feature maps');
set(gca,'YTick',size(weightImage,1)/2);
set(gca,'YTickLabelRotation',90);
colorbar;
set(gca,'FontSize',24,'FontName','Arial');
end

