offset = 1;
close all;
map = 1;
a = 100;
load([param.dataFolder param.affSubfolder param.affMaps(map).name '.mat']);
load([param.dataFolder param.outputSubfolder '\09072012-618765\seg1-2.mat']);
load('C:\minicube\testSeg2.mat');
load(param.cmSource);
%load([param.dataFolder param.outputSubfolder param.affMaps(map).name '\MorphRecon.mat']);
%% Affinity Map
img = imagesc(affX(5-offset:end-5-offset,5-offset:end-5-offset,a-offset), [-1.7 1.7]);
colormap('gray');
axis equal;
axis off;
axis tight;
saveas(gcf, 'aff.png');

%% Raw Data
figure;
imagesc(raw(5:end-5,5:end-5,a), [50 205]);
colormap('gray');
axis equal;
axis off;
saveas(gcf, 'raw.png');
%% Segmentation
figure;
b = label2rgb(v{4,4}(5-offset:end-5-offset,5-offset:end-5-offset,a-offset), autoKLEE_colormap, 'k');
imagesc(b);
axis equal;
axis off;
saveas(gcf, 'seg.png');
