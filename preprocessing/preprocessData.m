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
    clear vessels;
end
clear i j lastOfKnossosCube maskedRaw numberValidInCube raw vessels thisSliceBbox;
toc;

% Approximate gradients by first (in this function) mean downsampling
filterSize = [64; 64; 29];
[rawMean, x, y, z] = approximateGradients(vesselsMasked, dataset.bbox, filterSize);
save('/gaba/scratch/mberning/rawMean.mat');

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

display('Gradient correction');
thisSliceBbox = dataset.bbox;
[X,Y,Z] = meshgrid(y,x,z);
for i=1:length(zCoords)
    % Determine bounding box and interpolation for this z-slice
    thisSliceBbox(3,:) = [zCoords{i}(1) zCoords{i}(end)];	
	xq = thisSliceBbox(1,1):thisSliceBbox(1,2);
    yq = thisSliceBbox(2,1):thisSliceBbox(2,2);
    zq = thisSliceBbox(3,1):thisSliceBbox(3,2);
    [Xq,Yq,Zq] = meshgrid(yq,xq,zq);
    % Read original data and detected vessel
    % Interpolate correction voxel for each voxel and multiply
    correctionForSlice =  interp3(X,Y,Z,correctionVolume,Xq,Yq,Zq, 'linear', 121);
    warning off;
    raw = readKnossosRoi(vesselsMasked.root, vesselsMasked.prefix, thisSliceBbox); 
    vessels =  readKnossosRoi(strrep(vesselsMasked.root, '/color/', '/segmentation/'), vesselsMasked.prefix, thisSliceBbox, 'uint32', '', 'raw');
    warning on;
    raw = uint8(correctionForSlice .* double(raw));
	raw(vessels > 0) = 121;
    % Save to new datset to be used in pipeline
    warning off;
    writeKnossosRoi(gradientCorrected.root, gradientCorrected.prefix, thisSliceBbox(:,1)', raw);
    warning on;
    clear raw vessels; 
    Util.progressBar(i, length(zCoords));
end
clear X Xq Y Yq Z Zq;

% Create resoution pyramids for new dataset(s)
createResolutionPyramid(vesselsMasked.root);
createResolutionPyramid(gradientCorrected.root);
% Still a bit more complicated for downsampling segmentation(s)
thisRoot = strrep(vesselsMasked.root, '/color/', '/segmentation/');
createResolutionPyramid(thisRoot, vesselsMasked.prefix, dataset.bbox, strrep(thisRoot, '/1/', ''), true);

% Does not belong here (but who really cares, I do not anymore)
thisRoot = '/gaba/wKcubes/Connectomics department/sK15_Str_js_v3/segmentation/1/';
createResolutionPyramid(thisRoot, 'sK15_Str_js_v3_mag1', p.bbox, strrep(thisRoot, '/1/', ''), true);

