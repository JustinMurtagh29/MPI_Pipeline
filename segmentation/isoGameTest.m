% Can be used to generate uint8 stack of segmentation stored in variable
% seg
map = 5;
algo = 2;
r = 2;
par1 = 4;
par2 = 2;
load([param.dataFolder param.affSubfolder param.affMaps(map).name '.mat']);
load([param.dataFolder param.outputSubfolder param.affMaps(map).name '/seg' num2str(r) '-' num2str(algo) '.mat']);
seg = v{par1,par2};
clear v;

seg = single(seg);
sizeObj = hist(seg(:),1:max(seg(:)));
[sizeObjNew, order] = sort(sizeObj, 'descend');

% segCutoff = zeros(size(seg), 'single');
% for i=1:255
%     segCutoff(seg == order(i)) = i;
% end

% segCutoff = uint16(segCutoff);
% segNew = watershed_3D_PP( affX, affY, affZ, segCutoff);
aff = (affX + affY + affZ) ./ 3;
aff(aff < -0.8) = -0.8;
aff(aff > 1.2) = 1.2;
aff = aff - min(aff(:));
aff = aff ./ max(aff(:)).* 255;

%% Select all object belonging to a certain skeleton
load([param.dataFolder param.outputSubfolder param.affMaps(map).name '/evaluation' num2str(r) '-' num2str(algo) '.mat']);
splits = v.split(par1,par2).vec > 1;
idx = find(splits);

skeletonsToExtract = [30 52 53 68 74];
skelStack = cell(length(skeletonsToExtract),1);
for i=1:length(skeletonsToExtract)
    skelStack{i} = zeros(384,384,384);
    objects = v.split(par1,par2).obj{idx == skeletonsToExtract(i)};
    for j=1:length(objects)
        skelStack{i} = seg == objects(j) | skelStack{i};
    end
end

KLEE_v4('stack', raw, 'stack_2', uint8(skelStack{4}));

for i=1:length(skeletonsToExtract)
    skelStack{i} = 10 * bwdist(skelStack{i});
    writeKnossosCube(['/home/mberning/Desktop/isoTest/skel' num2str(i) '/'], 'segTest', [1 1 1], uint8(skelStack{i}(1:256,1:256,1:256)));
    system(['cp ' '/home/mberning/Desktop/isoTest/skel' num2str(i) '/x0001/y0001/z0001/segTest_x0001_y0001_z0001.raw ' '/home/mberning/Desktop/isoTest/skelWithDist' num2str(i) '.raw']);
end

%% Write data for transparency test

rawNorm = single(raw(1:256,1:256,1:256) - 100);
rawNorm(rawNorm < 0) = 0;
rawNorm(rawNorm > 80) = 80;
rawNorm = rawNorm / 80 * 255;
segOneObjID = seg(1:256,1:256,1:256);
segOneObjID(segOneObjID > 0) = 255;

KLEE_v4('stack', uint8(rawNorm), 'stack_2', uint8(segOneObjID));

writeKnossosCube('/home/mberning/Desktop/isoTest/temp/', 'segTest', [1 1 1], uint8(rawNorm));
system(['cp ' '/home/mberning/Desktop/isoTest/temp/x0001/y0001/z0001/segTest_x0001_y0001_z0001.raw ' '/home/mberning/Desktop/isoTest/rawNorm.raw']);
writeKnossosCube('/home/mberning/Desktop/isoTest/temp/', 'segTest', [1 1 1], uint8(segOneObjID));
system(['cp ' '/home/mberning/Desktop/isoTest/temp/x0001/y0001/z0001/segTest_x0001_y0001_z0001.raw ' '/home/mberning/Desktop/isoTest/segOneObjID.raw']);

%% Filter image with 3D-Gaussian blur
filter = fspecial3('gaussianFWHMhalfSize', 3);
rawBlurred = imfilter(rawNorm(1:256,1:256,1:256), filter, 'full');
KLEE_v4('stack', uint8(rawNorm), 'stack_2', uint8(raw));

%% writeKnossosRoi('/home/mberning/Desktop/isoTest/seg/', 'segTest', [1 1 1], uint8(segNew));
% writeKnossosRoi('/home/mberning/Desktop/isoTest/raw/', 'segTest', [1 1 1], raw);
raw512 = readKnossosRoi('/data/e_k0563/k0563_mag1/', '100527_k0563_mag1', [1000 1511; 1000 1511; 1000 1255]);
writeKnossosCube('/home/mberning/Desktop/isoTest/rawVeryBig/', 'segTest', [1 1 1], raw512);

%% writeKnossosRoi('/home/mberning/Desktop/isoTest/class/', 'segTest', [1 1 1], uint8(aff));
writeKnossosCube('/home/mberning/Desktop/isoTest/classBig/', 'segTest', [1 1 1], uint8(aff(1:256,1:256,1:256)));
writeKnossosCube('/home/mberning/Desktop/isoTest/rawBig/', 'segTest', [1 1 1], uint8(raw(1:256,1:256,1:256)));
addpath('KLEE');
% KLEE_v4('stack', raw, 'stack_2', seg, 'stack_3', segCutoff, 'stack_4', segNew);
KLEE_v4('stack', raw, 'stack_2', seg);
%% 
addpath('SVM');
border = uint8((seg == 0) * 255);
filt = fspecial3('gaussianFWHMhalfSize', 4);
borderBlurred = uint8(imfilter(single(border), filt));
borderExpanded = uint8((border ...
    | cat(1, zeros(1, size(border,2), size(border,3)), border(1:end-1,:,:)) ...
    | cat(1, border(2:end,:,:), zeros(1, size(border,2), size(border,3))) ...
    | cat(2, zeros(size(border,1), 1, size(border,3)), border(:,1:end-1,:)) ...
    | cat(2, border(:,2:end,:), zeros(size(border,1), 1, size(border,3))) ...
    | cat(3, border(:,:,2:end), zeros(size(border,1), size(border,2), 1)) ...
    | cat(3, zeros(size(border,1), size(border,2), 1), border(:,:,1:end-1)))*255);

mixed = raw;
mixed(borderBlurred > 0) = 50;

writeKnossosRoi('/home/mberning/Desktop/isoTest/border/', 'segTest', [1 1 1], border);
writeKnossosRoi('/home/mberning/Desktop/isoTest/borderBlurred/', 'segTest', [1 1 1], borderBlurred);
writeKnossosRoi('/home/mberning/Desktop/isoTest/borderExpanded/', 'segTest', [1 1 1], borderExpanded);
writeKnossosRoi('/home/mberning/Desktop/isoTest/mixed/', 'segTest', [1 1 1], mixed);
% 
addpath('KLEE');
KLEE_v4('stack', raw, 'stack_2', segNew, 'stack_3', border, 'stack_4', borderBlurred, 'stack_5', borderExpanded, 'stack_6', mixed );
% KLEE_v4('stack', raw(257:end,257:end,257:end), 'stack_2', borderBlurred(257:end,257:end,257:end));

