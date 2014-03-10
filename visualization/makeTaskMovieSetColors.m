function makeTaskMovieSetColors(segmentation, raw, outputFileAdress, colors)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

colormap = colors;

figure;
set(gcf,'NextPlot','replacechildren');
set(gcf,'Renderer','OpenGL');
writerObj = VideoWriter(outputFileAdress);
writerObj.FrameRate = 16;
open(writerObj);
for f=1:size(raw,3)
    hold off;
    imshow(raw(:,:,f), [60 180])
    hold on;
    temp = label2rgb(segmentation(:,:,f), colormap);
    himage = imshow(temp);
    set(himage, 'AlphaData', 0.2);
    drawnow;
    frame = getframe;
    writeVideo(writerObj,frame);
end
close(writerObj);
close all;

end