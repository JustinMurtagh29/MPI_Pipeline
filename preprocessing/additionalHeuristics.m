% Dataset to use for additional heuristics
dataset.saveFolder = '/gaba/scratch/mberning/temp/';
dataset.bbox_wK = [128, 128, 128, 5446, 8381, 3286];
dataset.experimentName = '2012-09-28_ex145_07x2_ROI2016_vessel';
dataset.raw.root = '/gaba/u/mberning/wkCubes/2012-09-28_ex145_07x2_ROI2016_corrected/color/1/';
dataset.raw.prefix = '2012-09-28_ex145_07x2_ROI2016_corrected_mag1';
dataset.raw.voxelSize = [11.24 11.24 28];
dataset = setParameterSettings(dataset);
dataset.seg.root = '/gaba/u/mberning/wkCubes/2012-09-28_ex145_07x2_ROI2016_vessel/segmentation/1/';
dataset.seg.prefix = '2012-09-28_ex145_07x2_ROI2016_vessel_mag1';

% Add border to each local bounding box
outsideBorder = [25 25 10];
desiredBorder = [200 200 80];
for x=1:numel(dataset.local)
    dataset.local(x).bboxBig = dataset.local(x).bboxSmall + [-desiredBorder; desiredBorder]';
    dataset.local(x).bboxBig(:,1) = max(dataset.bbox(:,1) - outsideBorder', dataset.local(x).bboxBig(:,1));
    dataset.local(x).bboxBig(:,2) = min(dataset.bbox(:,2) + outsideBorder', dataset.local(x).bboxBig(:,2));
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
    'cluster', '-l h_vmem=18G -l h_rt=24:00:00');
toc;
Cluster.waitForJob(job);

display('(Re)Downsampling KNOSSOS hierachies for segmentation to include updates');
tic;
thisBBox = [1 1 1; (ceil(dataset.bbox(:,2)./1024).*1024)']';
createResolutionPyramid(dataset.seg.root, dataset.seg.prefix, thisBBox, strrep(dataset.seg.root, '/1/', ''), true);
toc;

%% For debugging algorithm(s), look at results in webKnossos & add problematic locations here

% coord_wk = [1472, 5173, 220];
% coord_mat = coord_wk + 1;
 
% % Find linear indices where this data is processed and executed locally
% idx = cellfun(@(x)and(all(coord_mat >= x{1}{2}(:,1)'),all(coord_mat <= x{1}{2}(:,2)')), inputCell);
% theseInputs = cat(2, {dataset.raw, dataset.seg}, inputCell{idx}{1});
% dbstop in detectNucleiLocal at 11; 
% functionH(theseInputs{:});

