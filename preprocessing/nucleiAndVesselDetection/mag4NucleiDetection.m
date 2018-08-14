% Written by Christian Schramm <christian.schramm@brain.mpg.de>
% Script for applying nuclei detection in mag4 for a rough estimate to get
% the locations.

addpath(genpath('/u/mbeining/code/auxiliaryMethods/'));
addpath(genpath('/u/mbeining/code/pipeline/'));
outputFolder = '/tmpscratch/mbeining/data/cubing/MB_hc_l_31_st013_BKalign/nuclei/';

% Nuclei detection in mag4 to get their locations for following fine
% detection in mag1.
dataset.root = '/wKcubes_live/2018_08_06_MB_HCl31_st013_BKfull/color/1/';
dataset.backend = 'wkwrap';
dataset.blocksize = 32;
mag1bbox =  [1   28485;1   18320;1    512];%7164];

datasetMag4 = dataset;
datasetMag4.root = strrep(dataset.root,'/1/','/4-4-2/');

datasetMag4SegVessel = datasetMag4;
datasetMag4SegVessel.root = strrep(datasetMag4.root,'/color/','/segmentation/vessel/');
Util.mkdir(datasetMag4SegVessel.root)
datasetMag4SegNuclei = datasetMag4;
datasetMag4SegNuclei.root = strrep(datasetMag4.root,'/color/','/segmentation/nuclei/');
Util.mkdir(datasetMag4SegNuclei.root)
%% Load raw and vessels
mag4bbox = (mag1bbox - 1) ./4 + 1;
mag4bbox(:,1) = ceil(mag4bbox(:,1));
mag4bbox(:,2) = floor(mag4bbox(:,2));
raw = loadRawData(datasetMag4, mag4bbox);
if exist(datasetMag4SegVessel.root,'dir')
    loadSegDataGlobal(datasetMag4SegVessel, mag4bbox(:,1)');
else    
    % Applies the vessel detection
    vessels = detectVessels(raw, false);
    display('Start saving vessels');
    wkwInit('new', datasetMag4SegVessel.root, 32, 32, 'uint32',1)
    saveSegDataGlobal(datasetMag4SegVessel, mag4bbox(:,1)', uint32(vessels));
end
%% Detect nuclei with visualization flag set to false
% Set flag to true for new dataset parameter tuning
% The mask is necessary because of the zeros in the data (due to
% destruction of staining).
nuclei = detectNuclei(raw, vessels, false);

display('Start writing');
wkwInit('new', datasetMag4SegNuclei.root, 32, 32, 'uint32',1)
saveSegDataGlobal(datasetMag4SegNuclei, mag4bbox(:,1)', uint32(nuclei));

% Each connected area gets its own value, so each nucleus its own ID
nucleiL = labelmatrix(bwconncomp(nuclei));
% Determines centroid and extension of each nucleus
rp = regionprops(permute(nucleiL,[2 1 3]));  % this could be wrong

display('Start saving coordinates');
Util.save(fullfile(outputFolder,'NucleiCoordinates.mat'),rp)

