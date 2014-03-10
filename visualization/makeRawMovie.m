function makeRawMovie( segmentation, raw)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

load('C:\code\segmentation\autoKLEE_colormap.mat');
autoKLEE_colormap = repmat(autoKLEE_colormap, 100, 1);

figure;
set(gcf,'NextPlot','replacechildren');
set(gcf,'Renderer','OpenGL');
writerObj = VideoWriter('C:\Users\mberning\Desktop\bigRawMovie.avi');
writerObj.FrameRate = 32;
open(writerObj);
for f=1:size(raw,3)
    hold off;
    imshow(raw(:,:,f), 'InitialMagnification', 50)
    %imshow(raw(:,:,f), [60 180], 'InitialMagnification', 50);
    hold on;
    temp = label2rgb(segmentation(:,:,f), autoKLEE_colormap);
    himage = imshow(temp, [60 180], 'InitialMagnification', 50);
    set(himage, 'AlphaData', 0.15 );
    frame = getframe;
    writeVideo(writerObj,frame');
end
close(writerObj);
close all;

end