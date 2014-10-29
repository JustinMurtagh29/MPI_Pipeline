%% Set enviroment & parameters
clc; clear all;
addpath('auxiliary');
addpath('segmentation');
if strcmp(computer,'PCWIN64')
    param.dataFolder =  'C:\data\minicube\';
else
    param.dataFolder =  '/data/cortex/20130407cube/';
end
param.affSubfolder = 'aff/';
param.outputSubfolder = 'seg/';
param.figureSubfolder = [param.outputSubfolder 'figures' filesep];

% Set parameter for scan
param.r = 0; % Radii for Morphological Reconstruction
param.algo{1} = 'v1'; % Regionalmin Segmentation
param.pR{1,1} = {.18:0.01:0.24 0:100:200};
param.algo{2} = 'v2'; % Threshold Segmentation
param.pR{1,2} =  {-0.04:0.01:0.02 0:100:200}; % Threshold & marker volumes cutoff (should be increasing values for 2nd param)
% Set parameter for evaluation
param.nodeThres = 1; % Number of nodes within object that count as connection
param.skeletons = 'L4_dense_skeletons.726.nml'; % Skeleton file for segmentation evaluation
param.cmSource = 'segmentation/autoKLEE_colormap.mat';
param.makeSegMovies = true;
param.makeErrorMovies = true;
param.plotObjSizeHist = true;
param.objSizeCutoff = 100000;
param.plotObjChains = true;
param.plotSynRate = true;
param.makeErrorStacks = true;

% Some dependent parameter calculations (usually no change needed)
param.affMaps = dir([param.dataFolder param.affSubfolder '*.mat']);
for i=1:length(param.affMaps)
    param.affMaps(i).name(end-3:end) = [];
end
clear i;

%% Transform to local coordiantes & correct to BBox (careful with nodes & nodesNumDataAll)
skel = readKnossosNml([param.dataFolder param.skeletons]);
skel = switchToLocalCoords( skel, [4000 4500 2200], 0);
skel = correctSkeletonsToBBox(skel, [800 800 200]);
writeKnossosNml([param.dataFolder param.skeletons 'local'], skel);
param.skeletons = [param.skeletons 'local'];
skel = readKnossosNml([param.dataFolder param.skeletons]);
param.skel = skel;
param.totalPathLength = getPathLength(skel);

%% Save parameter file to disk
if ~exist([param.dataFolder param.outputSubfolder], 'dir')
    mkdir([param.dataFolder param.outputSubfolder]);
end
save([param.dataFolder param.outputSubfolder '/parameter.mat'], 'param');

%% Intermediate start
addpath('auxiliary');
addpath('segmentation');
subfolder = 'seg';
if strcmp(computer,'PCWIN64')
    dataFolder =  'C:\data\minicube\';
else
    dataFolder =  '/data/cortex/minicube/';
end
load([dataFolder subfolder '/parameter.mat']);

%% Main cell for parameter tuning
matlabpool 4;
% Perform morphological reconstruction (also complements; for r=0 just imcomplement is performed)
morphRcortex(param);
% Start the parameter scan (will automatically save & overwrite save to outputSubfolder)
scanParameterSegmentationCortex(param);
% Analysis
evalParameterSegmentation(param);
% Overview of performance of different segmentations
visualizeOverview_2(param);
matlabpool close;
