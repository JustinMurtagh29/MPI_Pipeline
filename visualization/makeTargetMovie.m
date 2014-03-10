function makeTargetMovie(target, raw, outputFileAdress)

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
    temp = label2rgb(target(:,:,f) ~= 1, [1 0 0]);
    himage = imshow(temp);
    set(himage, 'AlphaData', .2 );
    drawnow;
    frame = getframe;
    writeVideo(writerObj,frame);
end
close(writerObj);
close all;

end