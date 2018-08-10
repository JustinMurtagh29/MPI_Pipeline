% Settings
dataset.root = '/wKcubes_live/2018_08_06_MB_HCl31_st013_BKfull/color/1/';
dataset.backend = 'wkwrap';
dataset.bbox =  [1   28485;1   18320;1    512];%7164];
dataset.blocksize = 32;
% dataset.bbox = dataset.bbox + [-25 25; -25 25; -10 10];
datasetVesselsMasked.root = '/wKcubes_live/2018_08_06_MB_HCl31_st013_BKfull_vessel/color/1/';
datasetVesselsMasked.backend = 'wkwrap';
datasetVesselsMasked.blocksize = dataset.blocksize;
datasetGradientCorrected.root = '/wKcubes_live/2018_08_06_MB_HCl31_st013_corrected/color/1/';
datasetGradientCorrected.backend = 'wkwrap';
datasetGradientCorrected.blocksize = dataset.blocksize;

mkdir(datasetVesselsMasked.root)
mkdir(datasetGradientCorrected.root)
% Define slice to be loaded in memory in paralell (one KNOSSOS cube z-plane)
zCoords = dataset.bbox(3,1):dataset.bbox(3,2);
blocksizeToUse = dataset.blocksize *1;
lastOfWKWCube = [find(mod(zCoords,blocksizeToUse) == 0) length(zCoords)];
numberValidInCube = [lastOfWKWCube(1) diff(lastOfWKWCube)];
zCoords = mat2cell(zCoords, 1, numberValidInCube);
zCoords(numberValidInCube == 0) = [];

% Define regions (CC according to regionprops) to be deleted
% Done manually by looking at whether FP or TP vessel detection
% After some rounds I only checked > 10k voxel regions, all else seems wrong & unimportant
% To do so run once with all regions set to [], look at webKnossos at locations output in next for loop
% Errors were: Myelin, nuclei charging artifacts or close to dataset border

% Uncomment to get all detected vessel, disables selection above
BVregionsToRemove = cell(numel(zCoords),1);

% Detect vessels in each slice seperately, store raw data masked with mean at vessel locations and segmentation for detection visualization
display('Detecting blood vessels');
inputCell = cellfun(@(x,y) {[dataset.bbox(1:2,:);x(1) x(end)],y},zCoords,BVregionsToRemove','uni',0);

job = Cluster.startJob( ...
    @detectVesselsBBox, inputCell,'sharedInputs',{dataset,datasetVesselsMasked}, ...
    'cluster', {'memory', 32,'time', '50:00:00','taskConcurrency',20}, ...
    'name', 'vesselDetection','numOutputs',1);
Cluster.waitForJob(job);


display('Downsampling KNOSSOS hierachies');
tic;
% Create resoution pyramids for new dataset(s)
createResolutionPyramid(datasetVesselsMasked.root);
% Still a bit more complicated for downsampling segmentation(s)
thisRoot = strrep(datasetVesselsMasked.root, '/color/', '/segmentation/');
thisBBox = [1 1 1; (ceil(dataset.bbox(:,2)./1024).*1024)']';
createResolutionPyramid(thisRoot, datasetVesselsMasked.prefix, thisBBox, strrep(thisRoot, '/1/', ''), true);
toc;

% Approximate gradients by first (in this function) mean downsampling
filterSize = [64; 64; 29];
[rawMean, x, y, z] = approximateGradients(datasetVesselsMasked, dataset.bbox, filterSize);
Util.save('/gaba/scratch/mberning/rawMean.mat');

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
lastOfWKWCube = [find(mod(zCoords,32) == 0) length(zCoords)];
numberValidInCube = [lastOfWKWCube(1) diff(lastOfWKWCube)];
zCoords = mat2cell(zCoords, 1, numberValidInCube);
zCoords(numberValidInCube == 0) = [];

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
    raw = loadRawData(datasetVesselsMasked, thisSliceBbox); 
    vessels =  loadSegDataGlobal(struct( ...
        'root', strrep(datasetVesselsMasked.root, '/color/', '/segmentation/'), ...
        'prefix', datasetVesselsMasked.prefix), thisSliceBbox);
    raw = uint8(correctionForSlice .* double(raw));
	raw(vessels > 0) = 121;
    % Save to new datset to be used in pipeline
    saveRawData(datasetGradientCorrected, thisSliceBbox(:, 1)', raw);
    clear raw vessels;
    Util.progressBar(i, length(zCoords));
end
clear X Xq Y Yq Z Zq;

display('Downsampling KNOSSOS hierachies');
tic;
% Create resoution pyramids for new dataset(s)
createResolutionPyramid(datasetGradientCorrected.root);
toc;

