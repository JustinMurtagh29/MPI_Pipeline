function makeTaskMovieSet(raw, mito, outputFileAdress)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

figure;
set(gcf,'NextPlot','replacechildren');
set(gcf,'Renderer','OpenGL');
writerObj = VideoWriter(outputFileAdress);
writerObj.FrameRate = 16;
open(writerObj);
for f=1:size(raw,3)
    hold off;
    imshow(raw(:,:,f), [60 180]);
        hold on;
    temp = label2rgb(mito(:,:,f), [1 0 0]);
    himage = imshow(temp);
    set(himage, 'AlphaData', 0.1 );
    drawnow;
    frame = getframe;
    writeVideo(writerObj,frame);
end
close(writerObj);
close all;

end