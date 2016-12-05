% Dataset with raw data to use for additional heuristics
dataset.root = '/gaba/u/mberning/wkCubes/2012-09-28_ex145_07x2_ROI2016_corrected/color/1/';
dataset.prefix = '2012-09-28_ex145_07x2_ROI2016_corrected_mag1';
% Dataset with vessel segmentation (see preprocessData.m)
vesselsMasked.root = '/gaba/u/mberning/wkCubes/2012-09-28_ex145_07x2_ROI2016_vessel/segmentation/1/';
vesselsMasked.prefix = '2012-09-28_ex145_07x2_ROI2016_vessel_mag1';
% Define bounding box for detection
dataset.bbox = [128, 128, 128, 5573, 8508, 3413];
% "Training region nuclei" (where heuristics where optimized)
dataset.bbox = [2105, 2606, 740, 1000, 1000, 500];
% "Training region myelin" (where heuristics where optimized)
dataset.bbox = [128, 128, 128, 1000, 1000, 500];

% Due to new bounding box format
dataset.bbox(:,4:6) = dataset.bbox(:,1:3) + dataset.bbox(:,4:6) - 1;
dataset.bbox = Util.convertWebknossosToMatlabBbox(dataset.bbox);

% Define slice to be loaded in memory in paralell (one KNOSSOS cube z-plane)
zCoords = dataset.bbox(3,1):dataset.bbox(3,2);
lastOfKnossosCube = [find(mod(zCoords,128) == 0) length(zCoords)];
numberValidInCube = [lastOfKnossosCube(1) diff(lastOfKnossosCube)];
zCoords = mat2cell(zCoords, 1, numberValidInCube);
zCoords(numberValidInCube == 0) = [];

% Detect nuclei in (slightly overlapping) small virtual slices along z
borderToLoad = [25 25 10];
display('Detecting nuclei & myelin');
tic;
for i=1:length(zCoords)
    thisSliceBbox = dataset.bbox;
    thisSliceBbox(3,:) = [zCoords{i}(1) zCoords{i}(end)];
    thisSliceBboxWithBorder = thisSliceBbox + [-borderToLoad; borderToLoad]';
    % Load raw data and vessel from new dataset
    raw = readKnossosRoi(dataset.root, dataset.prefix, thisSliceBboxWithBorder);
    vessels = readKnossosRoi(vesselsMasked.root, vesselsMasked.prefix, thisSliceBbox, 'uint32');
    % Apply new heuristics
    nuclei = detectNucleiZSlice(raw, borderToLoad);
    myelin = detectMyelinZSlice(raw, borderToLoad);
    mito = detectMitoZSlice(raw, borderToLoad);
    % Exlusion heuristics, e.g. vessel might be detected as nuclei, myelin
    % as mitochondria
    nuclei(vessels) = 0;
    mito(myelin) = 0;
    % Make sure no voxels are labeled double before pasting together
        assert(~any((vessel+nuclei+myelin+mito) == 2));

    vessels(nuclei) = 2;
    writeKnossosRoi(vesselsMasked.root, vesselsMasked.prefix, thisSliceBbox(:,1)', uint32(vessels), 'uint32', '', 'noRead');
    Util.progressBar(i, length(zCoords));
    clear raw vessels nuclei;
end
clear i lastOfKnossosCube numberValidInCube raw vessels nuclei thisSliceBbox;
toc;

display('(Re)Downsampling KNOSSOS hierachies for segmentation to include updates (Nuclei, Myelin & Mitochondria)');
tic;
% Still a bit more complicated for downsampling segmentation(s)
thisRoot = strrep(vesselsMasked.root, '/color/', '/segmentation/');
thisBBox = [1 1 1; (ceil(dataset.bbox(:,2)./1024).*1024)']';
createResolutionPyramid(thisRoot, vesselsMasked.prefix, thisBBox, strrep(thisRoot, '/1/', ''), true);
toc;
