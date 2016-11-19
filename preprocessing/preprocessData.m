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
zCoords(numberValidInCube == 0) = [];

% Define regions (CC according to regionprops) to be deleted
% Done manually by looking at whether FP or TP vessel detection
% After some rounds I only checked > 10k voxel regions, all else seems wrong & unimportant
% To do so run once with all regions set to [], look at webKnossos at locations output in next for loop
% Errors were: Myelin, nuclei charging artifacts or close to dataset border
regions{1} = 4;
regions{2} = 4;
regions{4} = 5;
regions{7} = 6;
regions{8} = 6;
regions{9} = 6;
regions{10} = [5:9 11:13];
regions{13} = 7;
regions{15} = 4:6;
regions{16} = 3:6;
regions{17} = 3;
regions{18} = 4:13;
regions{19} = 3:33;
regions{20} = [1 4:15 17];
regions{21} = [4:10 12:20];
regions{22} = 4:44;
regions{23} = [2:4 7:40];
regions{24} = [3 6:30 32:73];
regions{26} = [5 6];
regions{27} = 4:8;
% Uncomment to get all detected vessel, disables selection above
%regions = cell(27,1);

% Detect vessels in each slice seperately, store raw data masked with mean at vessel locations and segmentation for detection visualization
thisSliceBbox = dataset.bbox;
display('Detecting blood vessels');
tic;
for i=1:length(zCoords) 
    thisSliceBbox(3,:) = [zCoords{i}(1) zCoords{i}(end)];
    raw = readKnossosRoi(dataset.root, dataset.prefix, thisSliceBbox);
    vessels = false(size(raw));
    for j=1:size(raw,3)
        vessels(:,:,j) = detectVesselsSingleImage(raw(:,:,j));
    end 
    % Regionprops CC of pixel (small compared to blood vessel)
    rp = regionprops(vessels, {'Area' 'Centroid' 'PixelIdxList'});
    display(['Iteration: ' num2str(i)]);
    for j=1:length(rp)
       display(['Region ' num2str(j) ': ' num2str(rp(j).Area) ' voxel, centroid: ' num2str(round(rp(j).Centroid([2 1 3]))+ thisSliceBbox(:,1)')]);
    end
    display('Regions removed:');
    display(num2str(regions{i}));
    vessels(cat(1,rp(regions{i}).PixelIdxList)) = 0;
    vessels = imclose(vessels, ones(1, 1, 7));  
    maskedRaw = raw;
    maskedRaw(vessels) = uint8(121);
    writeKnossosRoi(vesselsMasked.root, vesselsMasked.prefix, thisSliceBbox(:,1)', maskedRaw, 'uint8', '', 'noRead');
    writeKnossosRoi(strrep(vesselsMasked.root, '/color/', '/segmentation/'), vesselsMasked.prefix, thisSliceBbox(:,1)', uint32(vessels), 'uint32', '', 'noRead');
    %Util.progressBar(i, length(zCoords));
    clear vessels;
end
clear i j lastOfKnossosCube maskedRaw numberValidInCube raw vessels thisSliceBbox;
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
    raw = readKnossosRoi(vesselsMasked.root, vesselsMasked.prefix, thisSliceBbox); 
    vessels =  readKnossosRoi(strrep(vesselsMasked.root, '/color/', '/segmentation/'), vesselsMasked.prefix, thisSliceBbox, 'uint32');
    raw = uint8(correctionForSlice .* double(raw));
	raw(vessels > 0) = 121;
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

