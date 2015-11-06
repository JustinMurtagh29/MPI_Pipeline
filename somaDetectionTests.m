addpath(genpath('/home/mberning/code/manuelCode'));

root = '/home/mberning/fsHest/Data/berningm/20150205paper1submission/datasets/2012-09-28_ex145_07x2_corrected/color/4/';
prefix = '2012-09-28_ex145_07x2_mag4';

raw = readKnossosRoi(root, prefix, [1 2400; 1 1600; 1 860]);
save('/home/mberning/localStorage/l4_raw_mag4_temp.mat', '-v7.3');

%% Detect vessels with visualization flag set to false
% Set flag to true for new dataset parameter tuning
vessels = detectVessels(raw, false);

%% Make movie and isosurface visualization of detected vessels
makeSegMovie(vessels, raw, '/home/mberning/Desktop/vessels.avi', 1);

figure;
temp = vessels(1:6:end,1:6:end,:);
isosurface(temp, .5);
daspect([28 28 11.24]);
saveas(gcf, '/home/mberning/Desktop/vessels.png');

%% Detect vessels with visualization flag set to false
% Set flag to true for new dataset parameter tuning
somata = detectSomata(raw, vessels, true);

%% Make movie and isosurface visualization of detected vessels
makeSegMovie(somata, raw, '/home/mberning/Desktop/somata.avi', 1);

figure;
temp = somata(1:6:end,1:6:end,:);
isosurface(temp, .5);
daspect([28 28 11.24]);
saveas(gcf, '/home/mberning/Desktop/somata.png');
