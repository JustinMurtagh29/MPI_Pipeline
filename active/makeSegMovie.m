function makeSegMovie( segmentation, raw, filename)

load('/zdata/manuel/code/segmentation/autoKLEE_colormap.mat');
autoKLEE_colormap = repmat(autoKLEE_colormap, 100, 1);

figure('NextPlot', 'replacechildren', 'Visible', 'off');
%set(gcf,'Renderer','OpenGL');
writerObj = VideoWriter(filename, 'Uncompressed AVI');
writerObj.FrameRate = 4;
open(writerObj);
for f=1:size(raw,3)
    hold off;
    imshow(raw(:,:,f), [60 180]);
    hold on;
    temp = label2rgb(segmentation(:,:,f), autoKLEE_colormap);
    himage = imshow(temp, [60 180]);
    set(himage, 'AlphaData', 0.15 );
    frame = getframe;
    writeVideo(writerObj,frame);
end
close(writerObj);
close all;

end
