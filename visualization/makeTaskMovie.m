function makeTaskMovie( segmentation, raw, outputFileAdress)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

colormap = [1 0 0; 0 1 0];

figure;
set(gcf,'NextPlot','replacechildren');
set(gcf,'Renderer','OpenGL');
writerObj = VideoWriter(outputFileAdress);
writerObj.FrameRate = 16;
open(writerObj);
for f=1:size(raw,3)
    hold off;
    imshow(raw(:,:,f))
    hold on;
    temp = label2rgb(segmentation(:,:,f), colormap);
    himage = imshow(temp, [60 180]);
    set(himage, 'AlphaData', 0.15 );
    drawnow;
    frame = getframe;
    writeVideo(writerObj,frame);
end
close(writerObj);
close all;

end