% Load data
clear all; clc;
load I:\CortexConnectomics\Manuel\backup\20131126Laptop\Desktop\20130625allDataForFS.mat;

%%
x = 257:769;
y = 257:769;
z = 300;
aff = imcomplement(class);
fgm1 = imextendedmin(aff, .39, 26);
fgm2 = bwareaopen(fgm1, 10);

figure;
img = seg(x,y,z);
uniqueElem = unique(img);
for i=1:length(uniqueElem)
    img(img == uniqueElem(i)) = i - 1;
end
a = distinguishable_colors(length(unique(img)), [0 0 0; .1 .1 .1]);
idx = randperm(length(unique(img)));
a = a(idx,:);
img = label2rgb(img, a, [0 0 0]);
imagesc(img);
colormap(a);
axis equal; axis off;
saveas(gcf, 'C:\Users\mberning\Desktop\figureTemp\fig8h.tif');

figure;
img = raw(x,y,z);
imagesc(img);
caxis([40 200]);
colormap('gray');
axis equal; axis off;
saveas(gcf, 'C:\Users\mberning\Desktop\figureTemp\fig8c.tif');

figure;
img = class(x,y,z);
imagesc(img);
caxis([-1 1]);
colormap('gray');
axis equal; axis off;
saveas(gcf, 'C:\Users\mberning\Desktop\figureTemp\fig8d.tif');

%%
figure;
img = aff(x,y,z);
imagesc(img);
caxis([-0.5 2.5]);
colormap('gray');
axis equal; axis off;
saveas(gcf, 'C:\Users\mberning\Desktop\figureTemp\fig8e.tif');

figure;
img = seg(x,y,z);
for i=1:length(uniqueElem)
    img(img == uniqueElem(i)) = i - 1;
end
img(~fgm2(x,y,z)) = 0;
img = label2rgb(img, a, [0 0 0]);
imagesc(img);
colormap(a);
axis equal; axis off;
saveas(gcf, 'C:\Users\mberning\Desktop\figureTemp\fig8g.tif');

figure;
temp = label2rgb(fgm1(x,y,z) + fgm2(x,y,z), [1 0 0; 1 1 1], 'k');
imagesc(temp);
axis off;
axis equal;
saveas(gcf, 'C:\Users\mberning\Desktop\figureTemp\fig8f.png');
