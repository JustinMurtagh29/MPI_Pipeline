% Written by Christian Schramm <christian.schramm@brain.mpg.de>
% Script for applying nuclei detection in mag4 for a rough estimate to get
% the locations.

thresholdVessel = 200;
divideInN = 4;   % if dataset is too big, this number can be increased to chop the task into N parts
addpath(genpath('/gaba/u/mbeining/code/auxiliaryMethods/'));
addpath(genpath('/gaba/u/mbeining/code/pipeline/'));
outputFolder = '/tmpscratch/mbeining/data/cubing/MB_hc_l_31_st013_BKalign/nuclei/';
mkdir(outputFolder);
% Nuclei detection in mag4 to get their locations for following fine
% detection in mag1.
dataset.root = '/wKcubes_live/2018_08_06_MB_HCl31_st013_BKfull_adjHist/color/1/';
dataset.backend = 'wkwrap';
dataset.blocksize = 32;
mag1bbox =  [1   28485;1   18320;1    7164];
% mag1bbox =  [1   20000;1   10000;6401    6656];%7164];
datasetMag4 = dataset;
datasetMag4.root = strrep(dataset.root,'/1/','/4-4-2/');

datasetMag4Seg = datasetMag4;
datasetMag4Seg.root = strrep(datasetMag4.root,'/color/','/segmentation/heur/');
mkdir(datasetMag4Seg.root)


%% Load raw and vessels
mag4bbox = (mag1bbox - 1) ./4 + 1;
mag4bbox(:,1) = ceil(mag4bbox(:,1));
mag4bbox(:,2) = floor(mag4bbox(:,2));

zrange = floor(linspace(mag4bbox(3,1),mag4bbox(3,2),divideInN+1));
zrange = cell2mat(arrayfun(@(x,y) [x,y-1],zrange(1:end-1),zrange(2:end),'uni',0)');
zrange(end) = zrange(end) +1;
offset = zrange(1)-1;
part_mag4bbox = mag4bbox;
heurFull = uint8([]);
for s = 1:divideInN
    part_mag4bbox(3,:) = zrange(s,:);
    raw = loadRawData(datasetMag4, part_mag4bbox);
    try
        heur = loadSegDataGlobal(datasetMag4Seg, part_mag4bbox);
    end
    if ~exist('heur','var') || ~any(heur(:))
        % Applies the vessel detection
        display('Start vessel detection');
        heur = uint32(detectVessels(raw, thresholdVessel,4, 0));
        display('Start saving vessels');
        if ~exist(fullfile(datasetMag4Seg.root,'header.wkw'),'file')
            wkwInit('new', datasetMag4Seg.root, 32, 32, 'uint32',1)
        end
        saveSegDataGlobal(datasetMag4Seg, part_mag4bbox(:,1)', heur);
    end
    %% Detect nuclei with visualization flag set to false
    % Set flag to true for new dataset parameter tuning
    % The mask is necessary because of the zeros in the data (due to
    % destruction of staining).
    if ~any(heur(:) == 2)
        display('Start nuclei detection');
        heur = uint32(detectNuclei(raw, heur == 1, false,false)) * 2 + heur;
        
        display('Start writing');
        if ~exist(fullfile(datasetMag4Seg.root,'header.wkw'),'file')
            wkwInit('new', datasetMag4Seg.root, 32, 32, 'uint32',1)
        end
%         saveSegDataGlobal(datasetMag4Seg, part_mag4bbox(:,1)', heur);
    end
    heurFull = cat(3,heurFull,uint8(heur));
end
% finished parted heuristics, now load full mag4 stack into memory
clear raw heur
save('temp.mat','heurFull','-v7.3')
heurFull = bwareaopen(heurFull==2, 30000) * 2 + heurFull==1;
[x,y,z] = ind2sub(size(heurFull),find(heurFull == 2));
neighbors = bsxfun(@plus,[x,y,z],[-1,-1,-1;-1,-1,0;-1,0,0;0,-1,0;0,-1,-1;0,0,-1;-1,0,-1])
% stack is so big that uint32 would exceed memory of node, so do it in
% parts...
disp('saving nuclei heuristics')
for s = 1:divideInN
    part_mag4bbox(3,:) = zrange(s,:);
    saveSegDataGlobal(datasetMag4Seg, part_mag4bbox(:,1)', uint32(heurFull(:,:,zrange(s,1)-offset:zrange(s,2)-offset)));
end
% heur = loadSegDataGlobal(datasetMag4Seg, mag4bbox) == 2;
% Each connected area gets its own value, so each nucleus its own ID
cc = bwconncomp(heurFull == 2);
cc.PixelIdxList = cc.PixelIdxList(cellfun(@numel,cc.PixelIdxList)>30000);
cc.NumObjects = numel(cc.PixelIdxList);
nucleiL = labelmatrix(cc);
% Determines centroid and extension of each nucleus
rp = regionprops(permute(nucleiL,[2 1 3]));  % this could be wrong
display('Start saving coordinates');
Util.save(fullfile(outputFolder,'NucleiCoordinates.mat'),rp)