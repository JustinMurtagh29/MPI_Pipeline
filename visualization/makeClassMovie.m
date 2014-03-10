function makeClassMovie( class, raw )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

figure;
set(gcf,'NextPlot','replacechildren');
set(gcf,'Renderer','OpenGL');
writerObj = VideoWriter('C:\Users\mberning\Desktop\classMovie.avi');
writerObj.FrameRate = 16;
open(writerObj);
for f=1:size(raw,3)
    hold off;
    imagesc(raw(:,:,f));
    caxis([-2 2]);
    colormap('gray');
    hold on;
    temp = label2rgb(ones(size(class(:,:,f)),'uint8'), [1 0 0]);
    himage = imagesc(temp);
    set(gca, 'ALim', [0 1]);
    set(himage, 'AlphaDataMapping', 'scaled');
    set(himage, 'AlphaData', .2.*double(class(:,:,f)) );
    axis off;
    drawnow;
    frame = getframe;
    writeVideo(writerObj,frame);
    pause(.1);
end
close(writerObj);
close all;

end