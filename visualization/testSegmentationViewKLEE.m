%% Load raw data & affinity maps and generate segmentation
clear all;
clc;
load('/home/mberning/sync/fromFermat/toKLEE.mat');
rmpath('/home/mberning/code/KLEE/');
thres = affZ <= 1.71589994;
distTrans = bwdist(thres);

%% Watershed based (calculate forground markers first)
% for i=11:15
%     fgm{i} = imextendedmax(distTrans, i);
%     fgm{i} = imclose(fgm{i}, ones(5,5,5));
%     fgm{i} = imfill(fgm{i}, 'holes');
%     fgm{i} = bwareaopen(fgm{i}, 5);
%     distTransMod{i} = imimposemin(-distTrans, fgm{i});
%     segmentation2{i} = watershed(distTransMod{i},6);
% end

% fgm = imextendedmin(-distTrans, 15);
% fgm = imerode(fgm, ones(8,8,8));
% fgm = bwareaopen(fgm, 1000); 
% hminMinImposed = imimposemin(hmin, fgm);
rmpath('/home/mberning/code/KLEE/');
fgm = imcomplement(thres);
fgm2 = imfill(fgm,'holes');
fgm3 = imerode(fgm2, ones(9, 9, 9));
fgm4 = bwareaopen(fgm3, 50);
hmin = imhmin(-distTrans, .5);
hminMod = hmin;
hminMod(fgm3) = min(hmin(:));
test = watershed(hminMod);
addpath('/home/mberning/code/KLEE/');

%% View result in KLEE
addpath('/home/mberning/code/KLEE/');
KLEE_v4('stack', raw, 'stack_2', affX, 'stack_3', affY, 'stack_4', affZ);
KLEE_v4('stack', raw, 'stack_2', distTrans, 'stack_3', fgm3);
KLEE_v4('stack', raw, 'stack_2', test);
KLEE_v4('stack', hmin, 'stack_2', hminMod);
