addpath(genpath('/home/mberning/code/manuelCode'));

root = '/home/mberning/fsHest/Data/berningm/20150205paper1submission/datasets/2012-09-28_ex145_07x2_corrected/color/4/';
prefix = '2012-09-28_ex145_07x2_mag4';

raw = readKnossosRoi(root, prefix, [1 2400; 1 1600; 1 860]);
save('/home/mberning/localStorage/l4_raw_mag4_temp.mat', '-v7.3');

%% Detect vessels with visualization flag set to false
% Set flag to true for new dataset parameter tuning
vessels = detectVessels(raw, false);

%% Detect nuclei with visualization flag set to false
% Set flag to true for new dataset parameter tuning
nuclei = detectNuclei(raw, vessels, false);

%% Make movie and isosurface visualization of detected vessels
makeSegMovie(vessels, raw, '/home/mberning/Desktop/vessels.avi', 1);

figure;
temp = vessels(1:6:end,1:6:end,:);
isosurface(temp, .5);
daspect([28 28 11.24*6]);
saveas(gcf, '/home/mberning/Desktop/vessels.png');

%% Make movie and isosurface visualization of detected nuclei
makeSegMovie(nucleiPost, raw, '/home/mberning/Desktop/nuclei.avi', 1);

figure;
temp = nucleiPost(1:6:end,1:6:end,:);
isosurface(temp, .5);
daspect([28 28 11.24*6]);
saveas(gcf, '/home/mberning/Desktop/nuclei.png');

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

%% Write Knossos Hierachies

writeKnossosRoi(root, prefix, [1 2400; 1 1600; 1 860]);
