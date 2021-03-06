% Settings
addpath(genpath('/gaba/u/alik/code/pipeline'))
dataset.root = '/gaba/u/alik/wKcubes/2017-04-03_ppcAK99_76x79_ROI/color/1/';
dataset.prefix = '2017-04-03_ppcAK99_76x79_ROI_mag1';
dataset.bbox = [104,5983;104,7688;119,4822];%%full bbox ROI datasets=[1,6016; 1,7808;1,4864] WK=[0,0,0,6016,7808,4864]
%dataset.bbox = dataset.bbox + [-25 25; -25 25; -10 10];
vesselsMasked.root = '/gaba/u/alik/wKcubes/2017-04-03_ppcAK99_76x79_preprocessing/color/1/';
vesselsMasked.prefix = '2017-04-03_ppcAK99_76x79_preprocessing_mag1';
gradientCorrected.root = '/gaba/u/alik/wKcubes/2017-04-03_ppcAK99_76x79_corrected/color/1/';
gradientCorrected.prefix = '2017-04-03_ppcAK99_76x79_corrected_mag1';

% Define slice to be loaded in memory in paralell (one KNOSSOS cube z-plane)
zCoords = dataset.bbox(3,1):dataset.bbox(3,2);
lastOfKnossosCube = [find(mod(zCoords,128) == 0) length(zCoords)];
numberValidInCube = [lastOfKnossosCube(1) diff(lastOfKnossosCube)];
zCoords = mat2cell(zCoords, 1, numberValidInCube);
zCoords(numberValidInCube == 0) = [];

% Define regions (CC according to regionprops) to be deleted
% Done manually by looking at whether FP or TP vessel detection
% After some rounds I only checked > 10k voxel regions, all else seems wrong & unimportant
% To do so run once with all regions set to [], look at webKnossos at locations output in next for loop
% Errors were: Myelin, nuclei charging artifacts or close to dataset border
% regions{1} = 4;
% Uncomment to get all detected vessel, disables selection above
%regions = cell(27,1);

% Detect vessels in each slice seperately, store raw data masked with mean at vessel locations and segmentation for detection visualization
thisSliceBbox = dataset.bbox;
display('vesselMasking');
tic;
% Indices used for upsampling of labels to mag1 (start of if)
    x = ceil(0.25:0.25:((thisSliceBbox(1,2)-thisSliceBbox(1,1)+1)/4));
    x=[x(3:end),x(end),x(end)];
    y = ceil(0.25:0.25:((thisSliceBbox(2,2)-thisSliceBbox(2,1)+1)/4));
    y=[y(3:end),y(end),y(end)];
    
    
    
    load('/gaba/u/alik/scratch/cubing/2017-04-03_ppcAK99_76x79/results/vesselsFinal.mat');
for i=1:length(zCoords) 
    disp(num2str(i))
     thisSliceBbox(3,:) = [zCoords{i}(1) zCoords{i}(end)]
    
    z = ceil(0.25:0.25:((thisSliceBbox(3,2)-thisSliceBbox(3,1)+1)/4));
    % empirical correction
    z=[z(3:end),z(end),z(end)];
    
    raw = readKnossosRoi(dataset.root, dataset.prefix, thisSliceBbox);
    vesselsLocal= false(size(raw));
    thisBBoxMag4 = ceil(thisSliceBbox ./ 4)+ceil([1076-104, 1331-104,-119;1076-104, 1331-104,-119]./4)'
    theseVessels = vessels(thisBBoxMag4(1,1):thisBBoxMag4(1,2), ...
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
    % Regionprops CC of pixel (small compared to blood vessel)
    maskedRaw = raw;
    maskedRaw(vesselsLocal) = uint8(128);
    writeKnossosRoi(vesselsMasked.root, vesselsMasked.prefix,thisSliceBbox(:,1)', maskedRaw, 'uint8', '', 'noRead');
    writeKnossosRoi(strrep(vesselsMasked.root, '/color/', '/segmentation/'), vesselsMasked.prefix, thisSliceBbox(:,1)', uint32(vesselsLocal), 'uint32', '', 'noRead');
    %Util.progressBar(i, length(zCoords));
    
    clear vesselsLocal theseVessels z;
end
clear i j  maskedRaw raw  ;
toc;

display('Downsampling KNOSSOS hierachies');
tic;
% Create resoution pyramids for new dataset(s)
createResolutionPyramid(vesselsMasked.root);
% Still a bit more complicated for downsampling segmentation(s)
thisRoot = strrep(vesselsMasked.root, '/color/', '/segmentation/');
thisBBox = [1 1 1; (ceil(dataset.bbox(:,2)./1024).*1024)']';
createResolutionPyramid(thisRoot, vesselsMasked.prefix, thisBBox, strrep(thisRoot, '/1/', ''), true);
toc;

% Approximate gradients by first (in this function) mean downsampling
filterSize = [64; 64; 29];
[rawMean, x, y, z] = approximateGradients(vesselsMasked, dataset.bbox, filterSize);
Util.save('/gaba/scratch/alik/rawMean.mat');

% Gradient calculation
display('Smoothing gradient estimation');
tic;
% Calculate and smooth gradient along each axis (this will minimize effects
% of nuclei etc. variations on 10 micron scale as averaged over 100)
x_gradient = mean(mean(rawMean,2),3);
y_gradient = mean(mean(rawMean,1),3);
z_gradient = mean(mean(rawMean,1),2);
x_gradient_smooth = smooth(x_gradient, 5, 'lowess');
y_gradient_smooth = permute(smooth(y_gradient, 5, 'lowess'), [2 1 3]);
z_gradient_smooth = permute(smooth(z_gradient, 5, 'lowess'), [3 2 1]);
% Construct 3D volume with correction factors
correctionVolume = bsxfun(@times, x_gradient_smooth*y_gradient_smooth, z_gradient_smooth);
correctionVolume = mean(correctionVolume(:)) ./ correctionVolume;
% Extend correction volume (used for interpolation later) to avoid border
% effects within raw data
correctionVolume = padarray(correctionVolume, [2 2 2], 'symmetric');
dx = x(2)-x(1);
dy = y(2)-z(1);
dz = z(2)-z(1);
x = [x(1)-2*dx x(1)-dx x x(end)+dx x(end)+2*dx];
y = [y(1)-2*dy y(1)-dy y y(end)+dy y(end)+2*dy];
z = [z(1)-2*dz z(1)-dz z z(end)+dz z(end)+2*dz];
% Display some statistics to make sure everything is in a reasonable range
display(['Mean of correction: ' num2str(mean(correctionVolume(:)))]);
display(['Min of correction: ' num2str(min(correctionVolume(:)))]);
display(['Max of correction: ' num2str(max(correctionVolume(:)))]);
toc;

% Need to go to quarter KNOSSOS cube z-slices for next step, interpolation requires too much memory
zCoords = dataset.bbox(3,1):dataset.bbox(3,2);
lastOfKnossosCube = [find(mod(zCoords,32) == 0) length(zCoords)];
numberValidInCube = [lastOfKnossosCube(1) diff(lastOfKnossosCube)];
zCoords = mat2cell(zCoords, 1, numberValidInCube);
zCoords(numberValidInCube == 0) = [];

display('Gradient correction');
thisSliceBbox = dataset.bbox;
[X,Y,Z] = meshgrid(y,x,z);
for i=2:9
    % Determine bounding box and interpolation for this z-slice
    thisSliceBbox(3,:) = [zCoords{i}(1) zCoords{i}(end)];	
	xq = thisSliceBbox(1,1):thisSliceBbox(1,2);
    yq = thisSliceBbox(2,1):thisSliceBbox(2,2);
    zq = thisSliceBbox(3,1):thisSliceBbox(3,2);
    [Xq,Yq,Zq] = meshgrid(yq,xq,zq);
    % Read original data and detected vessel
    % Interpolate correction voxel for each voxel and multiply
    correctionForSlice =  interp3(X,Y,Z,correctionVolume,Xq,Yq,Zq, 'linear', 121);
    raw = readKnossosRoi(vesselsMasked.root, vesselsMasked.prefix, thisSliceBbox); 
    vessels =  readKnossosRoi(strrep(vesselsMasked.root, '/color/', '/segmentation/'), vesselsMasked.prefix, thisSliceBbox, 'uint32');
    raw = uint8(correctionForSlice .* double(raw));
	raw(vessels > 0) = 128;
    % Save to new datset to be used in pipeline
    writeKnossosRoi(gradientCorrected.root, gradientCorrected.prefix, thisSliceBbox(:,1)', raw);
    clear raw vessels;
    Util.progressBar(i, length(zCoords));
end
clear X Xq Y Yq Z Zq;

display('Downsampling KNOSSOS hierachies');
tic;
% Create resoution pyramids for new dataset(s)
createResolutionPyramid(gradientCorrected.root);
toc;
