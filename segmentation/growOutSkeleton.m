%% Parameter definition
addpath('auxiliary');
addpath('segmentation');
if strcmp(computer,'PCWIN64')
    param.dataFolder =  'C:\data\minicube\';
else
    param.dataFolder =  '/data/e_k0563/minicube/';
end
param.affSubfolder = ['affinityMaps' filesep];
param.outputSubfolder = ['growingOut' filesep];
param.figureSubfolder = [param.outputSubfolder 'figures' filesep];
param.skeletons = 'ssdf.203.nml';
param.cmSource = 'segmentation/autoKLEE_colormap.mat';
param.nodeThres = 1;
param.makeSegMovies = true;
param.makeErrorMovies = true;
param.plotObjSizeHist = true;
param.objSizeCutoff = 100000; % choose upper xlim bound histogram
param.plotObjChains = true;
param.plotSynRate = true;
param.makeErrorStacks = true;

%% Some dependent parameter calculations (usually no change needed)
param.affMaps = dir([param.dataFolder param.affSubfolder '*.mat']);
for i=1:length(param.affMaps)
    param.affMaps(i).name(end-3:end) = [];
end
clear i;

%% Some Preprocessing on the skeleton to keep KNOSSOS happy
skel = readKnossosNml([param.dataFolder param.skeletons]);
writeKnossosNml([param.dataFolder param.skeletons], skel);
skel = switchToLocalCoords( skel, [13 14 17] );
skel = correctSkeletonsToBBox(skel);
writeKnossosNml([param.dataFolder param.skeletons 'local'], skel);
param.skeletons = [param.skeletons 'local'];
skel = readKnossosNml([param.dataFolder param.skeletons]);
param.skel = skel;
param.totalPathLength = getPathLength(skel);

%% Transfer skeleton nodes to 3D matrix (for watershed) <-- sparse machen evtl sinnvoll
map = 4;
load([param.dataFolder param.affSubfolder param.affMaps(map).name '.mat']);
skelVol = zeros(size(affX));
for i=1:length(param.skel)
    for j=1:size(param.skel{i}.nodes, 1)
        skelVol(param.skel{i}.nodes(j,1), param.skel{i}.nodes(j,2), param.skel{i}.nodes(j,3)) = i;
    end
end
segMAT = imimposemin(affX, skelVol > 0);
segMAT = watershed(segMAT, 26);

map = 4;
load([param.dataFolder param.affSubfolder param.affMaps(map).name '.mat']);

seg = watershed_3D_PP(imcomplement(affX), imcomplement(affY), imcomplement(affZ), skelVol);

%% Anschauen in KLEE
addpath('KLEE');
KLEE_v4('stack', raw, 'stack_2', skelVol, 'stack_3', seg);
KLEE_v4('stack', raw, 'stack_2', imcomplement(affX), 'stack_3', imcomplement(affY), 'stack_4', imcomplement(affZ), 'stack_5', skelVol, 'stack_6', seg, 'stack_7', segMAT);
