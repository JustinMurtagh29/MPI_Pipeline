% Dataset to use for additional heuristics
dataset.saveFolder = '/gaba/scratch/mberning/temp/';
dataset.bbox_wK = [128, 128, 128, 5446, 8381, 3286];
dataset.experimentName = '2012-09-28_ex145_07x2_ROI2016_vessel';
dataset.raw.root = '/gaba/u/mberning/wkCubes/2012-09-28_ex145_07x2_ROI2016_vessel/color/1/';
dataset.raw.prefix = '2012-09-28_ex145_07x2_ROI2016_vessel_mag1';
dataset.raw.voxelSize = [11.24 11.24 28];
dataset = setParameterSettings(dataset);
dataset.seg.root = '/gaba/u/mberning/wkCubes/2012-09-28_ex145_07x2_ROI2016_vessel/segmentation/1/';
dataset.seg.prefix = '2012-09-28_ex145_07x2_ROI2016_vessel_mag1';

% Add fixed sized border to each local bounding box
borderToAdd = [25 25 10];
for x=1:numel(dataset.local)
    dataset.local(x).bboxBig = dataset.local(x).bboxSmall + [-borderToAdd; borderToAdd]';
end
dataset.seg = rmfield(dataset.seg, 'func');

% Detect nuclei & myelin in (overlapping) cubes
for i=1:numel(dataset.local)
    inputCell{i} = {{dataset.local(i).bboxBig, dataset.local(i).bboxSmall}};
end
functionH = @localDetection;
tic;
job = Cluster.startJob(functionH, inputCell, ...
    'name', 'nucleiDetection', ...
    'sharedInputs', {dataset.raw; dataset.seg}, ...
    'sharedInputsLocation', [1; 2], ...
    'cluster', '-l h_vmem=12G -l h_rt=24:00:00');
toc;
Cluster.waitForJob(job);

display('(Re)Downsampling KNOSSOS hierachies for segmentation to include updates');
tic;
thisBBox = [1 1 1; (ceil(dataset.bbox(:,2)./1024).*1024)']';
createResolutionPyramid(vesselsMasked.root, vesselsMasked.prefix, thisBBox, strrep(vesselsMasked.root, '/1/', ''), true);
toc;

