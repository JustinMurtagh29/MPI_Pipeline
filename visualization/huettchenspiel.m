addpath('C:\code\auxiliaryMethods\');
nx = 279;
ny = 126;
ddy = 65;
option = 'global';
y = 100;
z = 40;
speed = 0.1;
ddx = round(ddy/640*960);
dt = speed*ones(1,z);
dt(1:10) = linspace(.2, speed, 10);
% difficult1
%overallPos = [953 130 818];
% difficult2 
%overallPos = [1234, 414, 550];
% difficul3
%overallPos = [1521, 709, 222];
% difficul4
overallPos = [689, 475, 38]; 

%pos = [673 673+round(y/640*960); 293 293+y; 356 356+z];
pos = [overallPos(1)-round(y/640*960/2) overallPos(1)+round(y/640*960/2); ...
    overallPos(2)-round(y/2) overallPos(2)+round(y/2); ...
    overallPos(3)-round(z/2) overallPos(3)+round(z/2)];
a = double(readKnossosRoi('C:\2012-06-28_Cortex\mag1', '2012-06-28_Cortex_mag1', pos));
figure('Position', [520   378   560   420]); colormap('gray');
while true
    for i=1:z
        if strcmp(option, 'global') 
            imagesc(a(:,:,i)');
        else
            imagesc(a(nx-ddx:nx+ddx,ny-ddy:ny+ddy,i)');
        end
        axis off;
        if i==1
            pause;
        end
        %saveas(gcf, ['C:\Users\mberning\Desktop\movie\image' num2str(i, '%4i')], 'tif');
        imwrite(a(:,:,i)'/255, ['C:\Users\mberning\Desktop\movie\image' num2str(i, '%4i') '.tif'], 'tif');
        pause(dt(i));
        if i==z
            pause;
        end
    end
end