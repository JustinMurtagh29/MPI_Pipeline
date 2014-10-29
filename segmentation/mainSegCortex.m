%% Set enviroment & parameters
clc; clear param;
addpath(genpath('/data/segmentationOptimization/code/'));

% Folder structure 
param.dataFolder = '/data/segmentationOptimization/data/';
param.affSubfolder = 'aff/';
param.outputSubfolder = 'sixthRun/';
if ~exist([param.dataFolder param.outputSubfolder], 'dir');
    mkdir([param.dataFolder param.outputSubfolder])
end
param.figureSubfolder = [param.outputSubfolder 'figures' filesep];
if ~exist([param.dataFolder param.figureSubfolder], 'dir');
    mkdir([param.dataFolder param.figureSubfolder])
end

% Find all affinity maps from aff subfolder
param.affMaps = dir([param.dataFolder param.affSubfolder '/*.mat']);
for i=1:length(param.affMaps)
    param.affMaps(i).name(end-3:end) = [];
end
clear i;

% Set parameter for scan
param.r = 0; % Radii for Morphological Reconstruction
param.algo(1).fun = @(seg,pars) watershedSeg_v1_cortex(seg, pars(:));
param.algo(1).par = {0.02:0.02:0.7 0:50:100};
param.algo(2).fun = @(seg,pars) watershedSeg_v2_cortex(seg, pars(:));
param.algo(2).par = {};
param.algo(3).fun = @(seg,pars) watershedSeg_v3_cortex(seg, pars(:));
param.algo(3).par = {}; %if params are missing, algo is ignored
param.algo(4).fun = @(seg,pars) salient_watershedSeg_cortex(seg, pars(:));
param.algo(4).par = {}; %if params are missing, algo is ignored

% Set parameter for evaluation
param.nodeThres = 1; % Number of nodes within object that count as connection

% Set parameter for visualization of results
param.cmSource = 'seg/autoKLEE_colormap.mat';
param.makeSegMovies = true;
param.makeErrorMovies = true;
param.plotObjSizeHist = true;
param.objSizeCutoff = 100000;
param.plotObjChains = true;
param.plotSynRate = true;
param.makeErrorStacks = true;

%%  Choose affinity maps for training and test data
paramTest = param;
param.affMaps(2) = [];
paramTest.affMaps(1) = [];

%% Read skeleton for training
param.skeletons = 'L4_dense_skeletons.nml'; % Skeleton file for segmentation evaluation
param.skel = parseNml([param.dataFolder param.skeletons]);
param.skel = removeEmptySkeletons(param.skel);
% Switch to coordinates of small subcube (offset has to be inital voxel of
% bbox - [1 1 1] due to one index of matlab (another [1 1 1] if tracing was done in oxalis)
param.skel = switchToLocalCoords_v2(param.skel, [4097 4481 2250] - [1 1 1]);
% Remove all nodes outside small subcube
param.skel = correctSkeletonsToBBox_v2(param.skel, [640 768 201]);
% Calculate total path length of the skeleton within this box
param.totalPathLength = getPathLength(param.skel);
% Write skeleton video for control of training data
% load([param.dataFolder param.affSubfolder param.affMaps(1).name '.mat'], 'raw');
% raw = raw(1+25:end-25,1+25:end-25,1+10:end-10);
% makeSkeletonMovies(param, raw);
% Write local version of skeleton if needed
writeNml([param.dataFolder param.skeletons 'local'], param.skel);

%% Read skeleton for testing
paramTest.skeletons = 'DenseAlexcorrectedAnette.nml'; % Skeleton file for segmentation evaluation
paramTest.skel = parseNml([paramTest.dataFolder paramTest.skeletons]);
paramTest.skel = removeEmptySkeletons(paramTest.skel);
% Switch to coordinates of small subcube (offset has to be inital voxel of
% bbox - [1 1 1] due to one index of matlab (another [1 1 1] if tracing was done in oxalis)
paramTest.skel = switchToLocalCoords_v2(paramTest.skel, [1417 4739 890] - [1 1 1]);
% Remove all nodes outside small subcube
paramTest.skel = correctSkeletonsToBBox_v2(paramTest.skel, [300 300 300]);
% Calculate total path length of the skeleton within this box
paramTest.totalPathLength = getPathLength(paramTest.skel);
% Write skeleton video for control of testing data
% load([param.dataFolder param.affSubfolder param.affMaps(2).name '.mat'], 'raw');
% raw = raw(1+25:end-25,1+25:end-25,1+10:end-10);
% makeSkeletonMovies(param, raw);
writeNml([paramTest.dataFolder paramTest.skeletons 'local'], paramTest.skel);

%% Save parameter file to disk
save([param.dataFolder param.outputSubfolder 'parameter.mat'], 'param', 'paramTest');

%% Main cell for parameter tuning
matlabpool 2;
% Currently paralell over affinity maps, maybe change to parameter?
parfor map=1:length(param.affMaps)
    tic
    for r=1:length(param.r)
        disp(['Started morph, scan and eval for map ' num2str((map-1).*length(param.r)+r) '/' num2str(length(param.affMaps).*length(param.r))]);
        morphScanAndEval(param,param.affMaps(map).name,r);
    end
    toc
end
matlabpool close;

% Overview of performance of different segmentations
visualizeOverview_3(param);

%% Main cell for testing this shit
matlabpool 2;
% Currently paralell over affinity maps, maybe change to parameter?
parfor map=1:length(paramTest.affMaps)
    tic
    for r=1:length(paramTest.r)
        disp(['Started morph, scan and eval for map ' num2str((map-1).*length(paramTest.r)+r) '/' num2str(length(paramTest.affMaps).*length(paramTest.r))]);
        morphScanAndEval(paramTest,paramTest.affMaps(map).name,r);
    end
    toc
end
matlabpool close;

%% update 04.02
param.skel = parseNml([param.dataFolder param.skeletons]);
param.totalPathLength = getPathLength(param.skel);
param.skel = removeEmptySkeletons(param.skel);
paramTest.skel = parseNml([paramTest.dataFolder paramTest.skeletons]);
paramTest.totalPathLength = getPathLength(paramTest.skel);
paramTest.skel = removeEmptySkeletons(paramTest.skel);
save('/data/screen1_cortex.mat');
visualizeOverview_6cortex(param,paramTest);

%% update 25.02

load('/data/screen1_cortex.mat');
addpath(genpath('/data/segmentationOptimization/code/'));
% Execute cell 5 above to generate local version of the skeleton
% Then, edited equalizeSkeletons to just subsample (and write) skeleton 4
% (cortex test)
equalizeSkeletons();
paramTest.skeletons = [paramTest.skeletons 'local2'];
paramTest.skel = parseNml([paramTest.dataFolder paramTest.skeletons]);
paramTest.totalPathLength = getPathLength(paramTest.skel);
paramTest.skel = removeEmptySkeletons(paramTest.skel);
save('/data/screen1_cortexNew.mat');
visualizeOverview_6cortex(param,paramTest);
