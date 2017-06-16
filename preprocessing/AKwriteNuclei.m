addpath(genpath('/gaba/u/alik/code/pipeline'))
dataset.root = '/gaba/u/alik/wKcubes/2017-04-03_ppcAK99_76x79_ROI/color/1/';
dataset.prefix = '2017-04-03_ppcAK99_76x79_ROI_mag1';
dataset.bbox = [104,5983;104,7688;119,4822];
% Define slice to be loaded in memory in paralell (one KNOSSOS cube z-plane)
zCoords = dataset.bbox(3,1):dataset.bbox(3,2);
lastOfKnossosCube = [find(mod(zCoords,128) == 0) length(zCoords)];
numberValidInCube = [lastOfKnossosCube(1) diff(lastOfKnossosCube)];
zCoords = mat2cell(zCoords, 1, numberValidInCube);
zCoords(numberValidInCube == 0) = [];
if ~exist('nuclei','var')
    load('/gaba/scratch/alik/cubing/2017-04-03_ppcAK99_76x79/results/nucleiFinal.mat');
end
thisSliceBbox = dataset.bbox;
display('vesselMasking');
tic;
% Indices used for upsampling of labels to mag1 (start of if)
x = ceil(0.25:0.25:((thisSliceBbox(1,2)-thisSliceBbox(1,1)+1)/4));
x=[x(3:end),x(end),x(end)];
y = ceil(0.25:0.25:((thisSliceBbox(2,2)-thisSliceBbox(2,1)+1)/4));
y=[y(3:end),y(end),y(end)];


for i=1:length(zCoords)
    disp(num2str(i))
    thisSliceBbox(3,:) = [zCoords{i}(1) zCoords{i}(end)];
    
    z = ceil(0.25:0.25:((thisSliceBbox(3,2)-thisSliceBbox(3,1)+1)/4));
    % empirical correction
    z=[z(3:end),z(end),z(end)];
    
    
    thisBBoxMag4 = ceil(thisSliceBbox ./ 4)+ceil([1076-104, 1331-104,-119;1076-104, 1331-104,-119]./4)';
    theseVessels = nuclei(thisBBoxMag4(1,1):thisBBoxMag4(1,2), ...
        thisBBoxMag4(2,1):thisBBoxMag4(2,2), ...
        thisBBoxMag4(3,1):thisBBoxMag4(3,2));
    tic;
    if any(theseVessels(:))
        % Upsample labels to mag1
        vesselsLocal = theseVessels(x,y,z);
        se = strel('disk',5);
        vesselsLocal=imclose(vesselsLocal,se);
        
    end
    
    toc;

    
    writeKnossosRoi(strrep(dataset.root, '/color/', '/segmentation/'), dataset.prefix , thisSliceBbox(:,1)', uint32(vesselsLocal), 'uint32', '', 'noRead');
    %Util.progressBar(i, length(zCoords));
    
    clear vesselsLocal theseVessels z;
end

clear i j  maskedRaw raw  ;
toc;

display('Downsampling KNOSSOS hierachies');
tic;
thisRoot = strrep(dataset.root, '/color/', '/segmentation/');
thisBBox = [1 1 1; (ceil(dataset.bbox(:,2)./1024).*1024)']';
createResolutionPyramid(thisRoot, dataset.prefix, thisBBox, strrep(thisRoot, '/1/', ''), true);
toc;