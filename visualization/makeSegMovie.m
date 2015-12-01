function makeSegMovie( segmentation, raw, outputFileAdress, nrElementsInSegmentation)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4
    % -1 to not count zero as an element
    nrElementsInSegmentation = length(unique(segmentation(:))) - 1;
end

colors = distinguishable_colors(nrElementsInSegmentation, [0 0 0]);

figure;
set(gcf,'NextPlot','replacechildren');
set(gcf,'Renderer','OpenGL');
writerObj = VideoWriter(outputFileAdress);
writerObj.FrameRate = 16;
open(writerObj);
for f=1:size(raw,3)
    hold off;
    imshow(raw(:,:,f), [60 180], 'InitialMagnification', 50)
    hold on;
    temp = label2rgb(segmentation(:,:,f), colors, [1 1 1]);
    himage = imshow(temp);
    set(himage, 'AlphaData', 0.5);
    drawnow;
    frame = getframe;
    writeVideo(writerObj,frame);
end
close(writerObj);
close all;

end