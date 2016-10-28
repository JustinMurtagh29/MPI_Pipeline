% Settings
dataset.root = '/gaba/u/mberning/wkCubes/2012-09-28_ex145_07x2_ROI2016/color/1/';
dataset.prefix = '2012-09-28_ex145_07x2_ROI2016_mag1';
dataset.bbox = [128, 128, 128, 5573, 8508, 3413];
dataset.bbox = Util.convertWebknossosToMatlabBbox(dataset.bbox);
dataset.bbox = dataset.bbox + [-25 25; -25 25; -10 10];
vesselsMasked.root = '/gaba/u/mberning/wkCubes/2012-09-28_ex145_07x2_ROI2016_vessel/color/1/';
vesselsMasked.prefix = '2012-09-28_ex145_07x2_ROI2016_vessel_mag1';
gradientCorrected.root = '/gaba/u/mberning/wkCubes/2012-09-28_ex145_07x2_ROI2016_corrected/color/1/';
gradientCorrected.prefix = '2012-09-28_ex145_07x2_ROI2016_corrected_mag1';

% Define slice to be loaded in memory in paralell (one KNOSSOS cube z-plane)
zCoords = dataset.bbox(3,1):dataset.bbox(3,2);
lastOfKnossosCube = [find(mod(zCoords,128) == 0) length(zCoords)];
numberValidInCube = [lastOfKnossosCube(1) diff(lastOfKnossosCube)];
zCoords = mat2cell(zCoords, 1, numberValidInCube);

% Detect vessels in each slice seperately, store raw data masked with mean at vessel locations and segmentation for detection visualization
thisSliceBbox = dataset.bbox;
tic;
for i=1:length(zCoords) 
    thisSliceBbox(3,:) = [zCoords{i}(1) zCoords{i}(end)];
    warning off;
    raw = readKnossosRoi(dataset.root, dataset.prefix, thisSliceBbox);
    warning on;
    for j=1:size(raw,3)
        vessels(:,:,j) = detectVesselsSingleImage(raw(:,:,j));
    end 
    maskedRaw = raw;
    maskedRaw(vessels) = uint8(121);
    warning off;
    writeKnossosRoi(vesselsMasked.root, vesselsMasked.prefix, thisSliceBbox(:,1)', maskedRaw);
    writeKnossosRoi(strrep(vesselsMasked.root, '/color/', '/segmentation/'), vesselsMasked.prefix, thisSliceBbox(:,1)', uint32(vessels), 'uint32', '', 'noRead');
    warning on;
    Util.progressBar(i, length(zCoords));
end
clear i j;
toc;

% Calculate a mean downsampled, smoothed version of vesse
display('Gradient estimation');
tic;
filterSize = [64 64 29];
rawSmoothed = approximateGradients(vesselsMasked, dataset.bbox, filterSize);
toc;

display('Gradient correction');
tic;
for i=1:length(zCoords)
    thisSliceBbox(3,:) = [zCoords{i}(1) zCoords{i}(end)];
	z = ceil(i/filterSize(3));
	planeCorrectionImage = 1./(interp(mean(rawSmoothed(:,:,z),2),filterSize(1))*interp(mean(rawSmoothed(:,:,z),1),filterSize(2))./121^2);
	temp = readKnossosRoi(vesselsMasked.root, vesselsMasked.prefix, thisSliceBbox);
    vessels =  readKnossosRoi(strrep(vesselsMasked.root, '/color/', '/segmentation/'), vesselsMasked.prefix, thisSliceBbox, 'single', '', 'raw');
    temp = planeCorrectionImage.*temp;
	temp(vessels>0) = 121;
    writeKnossosRoi(gradientCorrected.root, gradientCorrected.prefix, thisSliceBbox(:,1)', temp);
    warning on
    Util.progressBar(length(zCoords), i);
end
toc;

