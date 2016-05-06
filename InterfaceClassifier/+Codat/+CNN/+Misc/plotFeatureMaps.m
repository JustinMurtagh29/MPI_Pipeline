function h = plotFeatureMaps( activity, layer, featureMaps )
%PLOTFEATUREMAPS Plot all feature maps in specified layer and input data.
% INPUT activity: The output of cnet.forward
%       layer: Integer specifying the layer.
%       featureMaps: Indices of feature maps as row vector. If not
%           specified then all feature maps of the corresponding layer are
%           plotted.
% OUTPUT h: Figure handle.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('featureMaps','var') || isempty(featureMaps)
    featureMaps = 1:size(activity{layer},4);
end
siz1 = size(activity{1},3);
sizL = size(activity{layer},3);
t = floor((siz1 - sizL)/2);
rows = ceil((length(featureMaps) + 1)/5);

%plot raw
subplot(rows,5,1);
imagesc(activity{1}(:,:,floor(sizL/2) + t,1));
j = 2;
colormap gray
axis off
title('Input data');

for fm = featureMaps
    subplot(rows,5,j);
    imagesc(activity{layer}(:,:,floor(sizL/2),fm));
    colormap gray
    axis off
    title(['Feature map ' num2str(fm)]);
    j = j + 1;
end
h = gcf;

end

