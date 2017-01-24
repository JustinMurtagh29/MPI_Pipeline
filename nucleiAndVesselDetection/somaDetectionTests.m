addpath(genpath('/home/mberning/code/manuelCode'));

root = '/home/mberning/fsHest/Data/berningm/20150205paper1submission/datasets/2012-09-28_ex145_07x2_corrected/color/4/';
prefix = '2012-09-28_ex145_07x2_mag4';

raw = readKnossosRoi(root, prefix, [1 2400; 1 1600; 1 860]);

%% Detect vessels with visualization flag set to false
% Set flag to true for new dataset parameter tuning
vessels = detectVessels(raw, false); % 7000s on 07x2 (mostly postprocessing)

%% Detect nuclei with visualization flag set to false
% Set flag to true for new dataset parameter tuning
nuclei = detectNuclei(raw, vessels, false);

%% Dilate with sphere to get outside nucleus and capture at least some endothelial segments
[x,y,z] = meshgrid(-3:3,-3:3,-1:1);
se = (x/3).^2 + (y/3).^2 + (z/1).^2 <= 1;
nuclei = imdilate(nuclei, se);
vessels = imdilate(vessels, se);

%% Make movie and isosurface visualization of detected vessels
makeSegMovie(vessels, raw, '/home/mberning/Desktop/vessels.avi', 1);

figure;
temp = vessels(1:6:end,1:6:end,:);
isosurface(temp, .5);
daspect([28 28 11.24*6]);
saveas(gcf, '/home/mberning/Desktop/vessels.png');

%% Make movie and isosurface visualization of detected nuclei
makeSegMovie(nuclei, raw, '/home/mberning/Desktop/nuclei.avi', 1);

figure;
temp = nuclei(1:6:end,1:6:end,:);
isosurface(temp, .5);
daspect([28 28 11.24*6]);
saveas(gcf, '/home/mberning/Desktop/nuclei.png');

%% Detect overlap
both = and(nuclei,vessels);

%% Make movie and isosurface visualization of detected overlaps
makeSegMovie(both, raw, '/home/mberning/Desktop/both.avi', 1);

figure;
temp = both(1:6:end,1:6:end,:);
isosurface(temp, .5);
daspect([28 28 11.24*6]);
saveas(gcf, '/home/mberning/Desktop/both.png');

%% Overlap should be vessels
vessels(both) = true;
nuclei(both) = false;

%% Save
clear x y z se both temp ans;
Util.save('/run/media/mberning/localStorage/07x2_mag4_vesselsAndNuclei.mat');

%% Joint isosurface

figure;
temp = nucleiPost(1:6:end,1:6:end,:);
n = patch(isosurface(temp, .5)); hold on;
temp = vessels(1:6:end,1:6:end,:);
v = patch(isosurface(temp, .5));

n.FaceColor = 'red';
v.FaceColor = 'blue';
n.EdgeColor = 'none';
v.EdgeColor = 'none';
daspect([28 28 11.24*6]);

view(3); axis tight
camlight 
lighting gouraud
saveas(gcf, '/home/mberning/Desktop/both.png');

%%
n = regionprops(nucleiPost);
v = regionprops(vessels);

figure;
hist([n(:).Area].*prod([28 28 11.24].*4)/1000^3, 50);
title('Nucleus size [microns^3]');

%%
rawNuclei = raw(nuclei);
regionprops(rawNuclei, 'all');

%% Give each nucleus its own unique ID
nucleiL = labelmatrix(bwconncomp(nuclei));
rp = regionprops(nucleiL);

%% Visualize all CC individually to judge/search errors
a = nucleiL(1:6:end,1:6:end,1:end);
for i=1:max(a(:))
    close all; figure;
    isosurface(a == i, .5);
    drawnow;
    saveas(gcf,['/home/mberning/Desktop/nuclei/' num2str(i, '%.3i') '_1.png']);
    camorbit(0,180); camlight;
    saveas(gcf,['/home/mberning/Desktop/nuclei/' num2str(i, '%.3i') '_2.png']);
    camorbit(90,0); camlight;
    saveas(gcf,['/home/mberning/Desktop/nuclei/' num2str(i, '%.3i') '_3.png']);
end

%% Visualize result with raw data for weird looking isosurfaces
% [9 68 77 96 111 117]
for i=117
    for z=floor(rp(i).BoundingBox(3)):ceil((rp(i).BoundingBox(3)+rp(i).BoundingBox(6)))
        x = floor(rp(i).BoundingBox(2)):ceil(rp(i).BoundingBox(2)+rp(i).BoundingBox(5));
        y = floor(rp(i).BoundingBox(1)):ceil(rp(i).BoundingBox(1)+rp(i).BoundingBox(4));
        imshow(raw(x,y,z)); hold on; 
        temp = label2rgb(nuclei(x,y,z), [0 1 0; 1 0 0; 0 0 1], [1 1 1]);
        himage = imshow(temp);
        set(himage, 'AlphaData', 0.3);
        title(num2str(i));
        pause(1);
        hold off;
    end
end

%% SAME POSTPROCESSING STARTS HERE, remove errors in detection

%% 9 is cell with two nuclei (glia) and one merged endothelial nucleus
i = 9;
x = floor(rp(i).BoundingBox(2)):ceil(rp(i).BoundingBox(2)+rp(i).BoundingBox(5));
y = floor(rp(i).BoundingBox(1)):ceil(rp(i).BoundingBox(1)+rp(i).BoundingBox(4));
z = floor(rp(i).BoundingBox(3)):ceil(rp(i).BoundingBox(3)+rp(i).BoundingBox(6));
isosurface(nuclei(x,y,z), .5);

% Test which size of structuring element we need to make it decompose into
% two parts
roi = nuclei(x,y,z);
distTraf = -bwdist(~roi);
minima = imextendedmin(distTraf, 1);
distTraf = imimposemin(distTraf, ~roi | minima);
fixed = watershed(distTraf);

% for i=1:size(fixed,3)
%     imagesc(fixed(:,:,i));
%     pause(.5);
% end

% Remove background object
fixed(fixed == 1) = 0;
% Seperate objects
roiNew = imclose(fixed == 2 | fixed == 3 | fixed == 4, ones(3,3,3)) | ...
    imclose(fixed == 5, ones(3,3,3));

% Fix in original data
nuclei(x,y,z) = roiNew;

%% 68 apical -> delete
nuclei(nucleiL == 68) = 0;

%% 111 merged blood vessel part
i = 111;
x = floor(rp(i).BoundingBox(2)):ceil(rp(i).BoundingBox(2)+rp(i).BoundingBox(5));
y = floor(rp(i).BoundingBox(1)):ceil(rp(i).BoundingBox(1)+rp(i).BoundingBox(4));
z = floor(rp(i).BoundingBox(3)):ceil(rp(i).BoundingBox(3)+rp(i).BoundingBox(6));

roi = nuclei(x,y,z);
distTraf = -bwdist(~roi);
minima = imextendedmin(distTraf, .7);
distTraf = imimposemin(distTraf, ~roi | minima);
fixed = watershed(distTraf);

for i=1:size(fixed,3)
    imagesc(fixed(:,:,i));
    title(num2str(i));
    pause(1);
end

fixed(fixed == 1) = 0;
fixed(fixed == 2) = 0;
roiNew = imclose(fixed == 3, ones(3,3,3));
nuclei(x,y,z) = roiNew;

%% 117 endothelial nucleus, but some apical attachted
i = 117;
x = floor(rp(i).BoundingBox(2)):ceil(rp(i).BoundingBox(2)+rp(i).BoundingBox(5));
y = floor(rp(i).BoundingBox(1)):ceil(rp(i).BoundingBox(1)+rp(i).BoundingBox(4));
z = floor(rp(i).BoundingBox(3)):ceil(rp(i).BoundingBox(3)+rp(i).BoundingBox(6));

roi = nuclei(x,y,z);
distTraf = -bwdist(~roi);
minima = imextendedmin(distTraf, 1);
distTraf = imimposemin(distTraf, ~roi | minima);
fixed = watershed(distTraf);

for i=1:size(fixed,3)
    imagesc(fixed(:,:,i));
    pause(.5);
end

fixed(fixed == 1) = 0;
fixed(fixed == 3) = 0;
roiNew = imclose(fixed == 2, ones(3,3,3));
nuclei(x,y,z) = roiNew;

%% Write Knossos Hierachies
nucleiL = labelmatrix(bwconncomp(nuclei));
rp = regionprops(nucleiL);

%%
a = nucleiL(1:6:end,1:6:end,1:end);
for i=1:max(a(:))
    close all; figure;
    isosurface(a == i, .5);
    drawnow;
    saveas(gcf,['/home/mberning/Desktop/nuclei/' num2str(i, '%.3i') '_1.png']);
    camorbit(0,180); camlight;
    saveas(gcf,['/home/mberning/Desktop/nuclei/' num2str(i, '%.3i') '_2.png']);
    camorbit(90,0); camlight;
    saveas(gcf,['/home/mberning/Desktop/nuclei/' num2str(i, '%.3i') '_3.png']);
end

%%
% [25 33 43 49 56 65 77 82 96 112 117 131 132]
for i=[65 77 82 96 112 117 131 132]
    for z=floor(rp(i).BoundingBox(3)):ceil((rp(i).BoundingBox(3)+rp(i).BoundingBox(6)))
        x = floor(rp(i).BoundingBox(2)):ceil(rp(i).BoundingBox(2)+rp(i).BoundingBox(5));
        y = floor(rp(i).BoundingBox(1)):ceil(rp(i).BoundingBox(1)+rp(i).BoundingBox(4));
        imshow(raw(x,y,z)); hold on; 
        temp = label2rgb(nuclei(x,y,z), [0 1 0; 1 0 0; 0 0 1], [1 1 1]);
        himage = imshow(temp);
        set(himage, 'AlphaData', 0.3);
        title(num2str(i));
        pause(.1);
        hold off;
    end
end

nuclei(nucleiL == 65) = 0;

%%
nucleiL = labelmatrix(bwconncomp(nuclei));
rp = regionprops(nucleiL);

%%
writeKnossosRoi('/run/media/mberning/localStorage/07x2_mag4_vessels/', 'vessels', [1 1 1], vessels);
writeKnossosRoi('/run/media/mberning/localStorage/07x2_mag4_nuclei/', 'nuclei', [1 1 1], nucleiL);
