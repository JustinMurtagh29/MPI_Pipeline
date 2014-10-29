% %% Set enviroment & parameters
% clc;
% addpath('auxiliary');
% addpath('segmentation');
% if strcmp(computer,'PCWIN64')
%     param.dataFolder =  'C:\data\minicube2extra\';
% else
%     param.dataFolder =  '/data/e_k0563/minicube2extra/';
% end
% param.affSubfolder = ['affinityMaps' filesep];
% param.outputSubfolder = ['seg' date filesep];
% param.figureSubfolder = [param.outputSubfolder 'figures' filesep];
% 
% % Set parameter for scan
% param.r = 0:1; % Radii for Morphological Reconstruction
% param.algo{1} = 'v1'; % Hmin Segmentation (not used for now)
% param.pR{1,1} = {[] []}; % Depth minima removal & marker volumes cutoff (should be increasing values for 2nd param)
% param.pR{2,1} = {[] []};
% param.algo{2} = 'v2'; % Threshold Segmentation
% param.pR{1,2} = {-0.16 20}; % Threshold & marker volumes cutoff (should be increasing values for 2nd param)
% param.pR{2,2} = {0.31 20};
% 
% % Set parameter for evaluation
% param.nodeThres = 1; % Number of nodes within object that count as connection
% param.skeletons = 'mini2eClean.nml'; % Skeleton file for segmentation evaluation
% param.cmSource = 'segmentation/autoKLEE_colormap.mat';
% param.makeSegMovies = true;
% param.makeErrorMovies = true;
% param.plotObjSizeHist = true;
% param.objSizeCutoff = 100000; % choose upper xlim bound histogram
% param.plotObjChains = true;
% param.plotSynRate = true;
% param.makeErrorStacks = true;
% 
% % Some dependent parameter calculations (usually no change needed)
% param.affMaps = dir([param.dataFolder param.affSubfolder '*.mat']);
% for i=1:length(param.affMaps)
%     param.affMaps(i).name(end-3:end) = [];
% end
% clear i;

% skel = readKnossosNml([param.dataFolder param.skeletons]);
% toDelete = []; % Collect empty skeletons
% for l=2:size(skel,2)
%     if ~isfield(skel{l}, 'nodesNumDataAll')
%         toDelete = [toDelete l];
%     end
% end
% parameters = skel{1}.parameters;
% skel(toDelete) = [];
% skel{1}.branchpoints = [];
% skel{1}.branchpointString = {};
% skel{1}.parameters = parameters;
% skel{1}.commentsString = {'<comments> </comments>'};
% writeKnossosNml([param.dataFolder param.skeletons 'cleanedUp'], skel);
% 
% %% Transform to local coordiantes & correct to BBox (careful with nodes & nodesNumDataAll)
% skel = readKnossosNml([param.dataFolder param.skeletons]);
% writeKnossosNml([param.dataFolder param.skeletons], skel);
% skel = switchToLocalCoords( skel, [24 4 4] );
% skel = correctSkeletonsToBBox(skel);
% writeKnossosNml([param.dataFolder param.skeletons 'local'], skel);
% param.skeletons = [param.skeletons 'local'];
% skel = readKnossosNml([param.dataFolder param.skeletons]);
% param.skel = skel;
% param.totalPathLength = getPathLength(skel);
% 
% %% Save parameter file to disk
% if ~exist([param.dataFolder param.outputSubfolder], 'dir')
%     mkdir([param.dataFolder param.outputSubfolder]);
% end
% save([param.dataFolder param.outputSubfolder 'parameter.mat'], 'param');

% Load old parameter file from disk
dateToLoad = 'seg14-Nov-2012';
addpath('auxiliary');
addpath('segmentation');
subfolder = ['seg' dateToLoad filesep];
if strcmp(computer,'PCWIN64')
    dataFolder =  'C:\data\minicube\';
else
    dataFolder =  '/data/e_k0563/minicube/';
end
load([dataFolder subfolder 'parameter.mat']);

% %% Main cell for parameter tuning
% matlabpool 3;
% % Perform morphological reconstruction (also complements; for r=0 just imcomplement is performed)
% morphR(param);
% % Start the parameter scan (will automatically save & overwrite save to outputSubfolder)
% scanParameterSegmentation(param);
% % Analysis
 evalParameterSegmentation(param);
% % Overview of performance of different segmentations
% visualizeOverview(param);
% matlabpool close;

%% Plot everything (including movies) for a specific segmentation
map = [1];
algo = [2];
r = [1];
par1 = [1];
par2 = [1];
for i=1:1
    visualizeSingle(param, map(i), algo(i), r(i), par1(i), par2(i));
end

map = [2];
algo = [2];
r = [2];
par1 = [1];
par2 = [1];
for i=1:1
    visualizeSingle(param, map(i), algo(i), r(i), par1(i), par2(i));
end

