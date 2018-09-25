% Written by Christian Schramm <christian.schramm@brain.mpg.de>
% and Marcel Beining <marcel.beining@brain.mpg.de>
% Script for applying vessel and  nuclei detection in mag4 for a rough estimate to get
% the nuclei locations.
% the heuristic mag4 stack contains the BVs (value 128) and the nuclei
% (value 255)
thresholdVessel = 200;
divideInN = 4;   % if dataset is too big, this number can be increased to chop the task into N parts
addpath(genpath('/gaba/u/mbeining/code/auxiliaryMethods/'));
addpath(genpath('/gaba/u/mbeining/code/pipeline/'));
outputFolder = '/tmpscratch/mbeining/data/cubing/MB_hc_l_31_st013_BKalign/nuclei/';
mkdir(outputFolder);
% Nuclei detection in mag4 to get their locations for following fine
% detection in mag1.
dataset.root = '/wKlive/2018_08_06_MB_HCl31_st013_BKfull_adjHist/color/1/';
dataset.backend = 'wkwrap';
dataset.blocksize = 32;
mag1bbox =  [1   28485;1   18320;1    7164];
% mag1bbox =  [1   20000;1   10000;6401    6656];%7164];
datasetMag4 = dataset;
datasetMag4.root = strrep(dataset.root,'/1/','/4-4-2/');


datasetMag4Heur = datasetMag4;
datasetMag4Heur.root = strrep(datasetMag4.root,'/color/','/heur/');
mkdir(datasetMag4Heur.root)

%% Load raw and vessels
mag4bbox = (mag1bbox - 1) ./4 + 1;
mag4bbox(:,1) = ceil(mag4bbox(:,1));
mag4bbox(:,2) = floor(mag4bbox(:,2));

zrange = floor(linspace(mag4bbox(3,1),mag4bbox(3,2),divideInN+1));
zrange = cell2mat(arrayfun(@(x,y) [x,y-1],zrange(1:end-1),zrange(2:end),'uni',0)');
zrange(end) = zrange(end) +1;
offset = zrange(1)-1;
part_mag4bbox = mag4bbox;

if ~exist(fullfile(datasetMag4Heur.root,'header.wkw'),'file')
    wkwInit('new', datasetMag4Heur.root, 32, 32, 'uint8',1)
    display('Start vessel detection');
    for s = 1:divideInN
        part_mag4bbox(3,:) = zrange(s,:);
        raw = loadRawData(datasetMag4, part_mag4bbox);
        % Applies the vessel detection
        heur = uint8(detectVessels(raw, thresholdVessel,4, 0))*128;
        %sanity check
        if any(all(all(heur,2),1))
            error('Some plane was fully labeled as blood vessel')
        end
        display('Start saving part of vessels');
        saveRawData(datasetMag4Heur, part_mag4bbox(:,1)', uint8(heur));
    end
end
display('Finished vessel detection');
cubesize = 800;
shift = 448;
% produce all cube starting points but make two non overlapping groups of
% them so that each can run in parallel without cube overwriting problems
[X,Y,Z] = meshgrid(mag4bbox(1,1):2*shift:mag4bbox(1,2),mag4bbox(2,1):2*shift:mag4bbox(2,2),mag4bbox(3,1):2*shift:mag4bbox(3,2));
offsets1 = mat2cell([X(:),Y(:),Z(:)],ones(numel(X),1),3);
[X,Y,Z] = meshgrid(mag4bbox(1,1)+shift:2*shift:mag4bbox(1,2),mag4bbox(2,1)+shift:2*shift:mag4bbox(2,2),mag4bbox(3,1)+shift:2*shift:mag4bbox(3,2));
offsets2 = mat2cell([X(:),Y(:),Z(:)],ones(numel(X),1),3);

display('Start nuclei detection');
job = Cluster.startJob( ...
    @detectNucleiInCube, offsets1,'sharedInputs',{datasetMag4, datasetMag4Heur,cubesize,true}, ...
    'cluster', {'memory', 12,'time', '30:00:00','taskConcurrency',20,'scheduler','slurm'}, ...
    'name', 'nucleiDetection');
Cluster.waitForJob(job);

job = Cluster.startJob( ...
    @detectNucleiInCube, offsets2,'sharedInputs',{datasetMag4, datasetMag4Heur,cubesize,true}, ...
    'cluster', {'memory', 12,'time', '30:00:00','taskConcurrency',20,'scheduler','slurm'}, ...
    'name', 'nucleiDetection');
Cluster.waitForJob(job);

heur = loadRawData(datasetMag4Heur, mag4bbox);
% raw = loadRawData(datasetMag4, mag4bbox);
%% Detect nuclei with visualization flag set to false
% Set flag to true for new dataset parameter tuning
% The mask is necessary because of the zeros in the data (due to
% destruction of staining).

% Each connected area gets its own value, so each nucleus its own ID
cc = bwconncomp(heur == 2);
cc.PixelIdxList = cc.PixelIdxList(cellfun(@numel,cc.PixelIdxList)>30000);
cc.NumObjects = numel(cc.PixelIdxList);
nucleiL = labelmatrix(cc);
% Determines centroid and extension of each nucleus
rp = regionprops(nucleiL);  
display('Start saving coordinates');
Util.save(fullfile(outputFolder,'NucleiCoordinates.mat'),rp)