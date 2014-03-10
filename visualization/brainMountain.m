addpath('C:\code\auxiliaryMethods\');
xy = 300;
z = 200;
dt = .1;
pos = [673 673+xy; 293 293+xy; 356 356+z];
a = readKnossosRoi('C:\2012-06-28_Cortex\mag1', '2012-06-28_Cortex_mag1', pos);
a = -double(a) + 255;
b = smooth3(double(a(:,:,:)), 'bbox', [5 5 5]);
for i=1:z+1
    b(:,:,i) = b(:,:,i) - mean(reshape(b(:,:,i),(xy+1)*(xy+1),1));
    b(:,:,i) = b(:,:,i) / std(reshape(b(:,:,i),(xy+1)*(xy+1),1));
end
close all;
writerObj = VideoWriter('C:\Users\mberning\Videos\mountains.avi');
writerObj.FrameRate = 4;
open(writerObj);
figure('Position', [1 41 1920 1084]);
cm = colormap('hot');
colormap(flipud(cm));
az = -55;
el = 78;
for i=1:z+1;
    h = surf(a(:,:,i));
    set(gca, 'CLim', [20 240]);
    set(h, 'edgecolor', 'none');
    axis off;
    view(az, el);
    shading interp;
    if i == 1
        pause;
        tic;
    end
    drawnow; 
    writeVideo(writerObj,getframe(gca, [1 1 1480 880]));
    pause(.1);
end;
close(writerObj);
