%% Load data

clear all; clc;
load /p/testSegControlRegion.mat;

%% Morphological Reconstruction

rT=1;
[x,y,z] = meshgrid(-rT:rT,-rT:rT,-rT:rT);
se = (x/rT).^2 + (y/rT).^2 + (z/rT).^2 <= 1;
% Opening by reconstruction
affEroded = imerode(aff, se);
affRecon = imreconstruct(affEroded, aff);
% Closing by reconstruction
affReconDilated = imdilate(affRecon, se);
affReconRecon= imreconstruct(imcomplement(affReconDilated), imcomplement(affRecon));

%% Marker

fgm1 = imregionalmin(affReconRecon);
fgm2 = bwareaopen(fgm1, 100, 26);

map1 = imimposemin(affReconRecon, fgm2);
seg1 = watershed(map1, 26);

bw1 = affReconRecon < 0.2;
bw1 = imopen(bw1, se);
bw2 = bwareaopen(bw1, 100, 26);

map2 = imimposemin(affReconRecon, bw2);
seg2 = watershed(map2, 26);

%% Visualization

% 1st row: Raw & Affinities
subplot(3,4,1);
imagesc(raw(:,:,100));
title('Raw');
colormap('gray');
axis off;
axis equal;
subplot(3,4,2);
imagesc(-aff(:,:,100));
title('Aff');
colormap('gray');
axis off;
axis equal;
subplot(3,4,3);
imagesc(-affRecon(:,:,100));
title('Opening by Reconstruction');
colormap('gray');
axis off;
axis equal;
subplot(3,4,4);
imagesc(affReconRecon(:,:,100));
title('Opening & Closing by Reconstruction');
colormap('gray');
axis off;
axis equal;
% Approach 1
subplot(3,4,6);
imagesc(fgm1(:,:,100));
title('Intermediate Foreground(regionalmin)');
colormap('gray');
axis off;
axis equal;
subplot(3,4,7);
imagesc(fgm2(:,:,100));
title('Foreground Marker 1(regionalmin + areaopen)');
colormap('gray');
axis off;
axis equal;
subplot(3,4,8);
imagesc(map1(:,:,100));
title('Input to Watershed 1');
colormap('gray');
axis off;
axis equal;
subplot(3,4,5);
imagesc(seg1(:,:,100));
title('Segmentation Approach 1');
colormap('gray');
axis off;
axis equal;
% Approach 2
subplot(3,4,10);
imagesc(bw1(:,:,100));
title('Intermediate Foreground (threshold + opening)');
colormap('gray');
axis off;
axis equal;
subplot(3,4,11);
imagesc(bw2(:,:,100));
title('Foreground Marker 2(threshold + areaopen)');
colormap('gray');
axis off;
axis equal;
subplot(3,4,12);
imagesc(map2(:,:,100));
title('Input to Watershed 2');
colormap('gray');
axis off;
axis equal;
subplot(3,4,9);
imagesc(seg2(:,:,100));
title('Segmentation Approach 2');
colormap('gray');
axis off;
axis equal;
