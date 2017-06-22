% Dataset to use for additional heuristics
dataset.saveFolder = '/gaba/scratch/alik/temp/';
dataset.bbox_wK = [128, 128, 128, 5850, 7535, 4684];%[103, 103, 118, 5880, 7585, 4704]
dataset.experimentName = '2017-04-03_ppcAK99_76x79_corrected';
dataset.raw.root = '/gaba/u/alik/wKcubes/2017-04-03_ppcAK99_76x79_corrected/color/1/';
dataset.raw.prefix = '2017-04-03_ppcAK99_76x79_corrected_mag1';
dataset.raw.voxelSize = [12 12 30];
dataset = setParameterSettings(dataset);
dataset.seg.root = '/gaba/u/alik/wKcubes/2017-04-03_ppcAK99_76x79_preprocessing/segmentation/1/';
dataset.seg.prefix = '2017-04-03_ppcAK99_76x79_preprocessing_mag1';

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
    inputCell{i} = {dataset.local(i).bboxBig, dataset.local(i).bboxSmall};
end
% Chose whether to run nuclei or myelin detection alone
functionH = @localDetectionMyelin;
%functionH = @localDetectionNuclei;

%% Job submissioniii
tic;
job = Cluster.startJob(functionH, inputCell, ...
    'name', 'myelinDetection', ...
    'sharedInputs', {dataset.raw; dataset.seg}, ...
    'sharedInputsLocation', [1; 2], ...
    'cluster', '-l h_vmem=18G -l h_rt=24:00:00 -tc 40');
toc;
Cluster.waitForJob(job);

display('(Re)Downsampling KNOSSOS hierachies for segmentation to include updates');
tic;
    gradientCorrected.root = '/gaba/u/alik/wKcubes/2017-04-03_ppcAK99_76x79_corrected/segmentation/1/';
    gradientCorrected.prefix = '2017-04-03_ppcAK99_76x79_corrected_mag1';
thisBBox = [1 1 1; (ceil(dataset.bbox(:,2)./1024).*1024)']';
createResolutionPyramid(gradientCorrected.root, gradientCorrected.prefix, thisBBox, strrep(gradientCorrected.root, '/1/', ''), true);
toc;
%% For debugging algorithm(s), look at results in webKnossos & add problematic locations here
% 
% coord_wk = [3008, 6970, 2902];
% coord_mat = coord_wk + 1;
% 
% % Find linear indices where this data is processed and executed locally
% idx = cellfun(@(x)and(all(coord_mat >= x{1}{2}(:,1)'),all(coord_mat <= x{1}{2}(:,2)')), inputCell);
% theseInputs = cat(2, {dataset.raw, dataset.seg}, inputCell{idx}{1});
% 
% functionH(theseInputs{:});
% % Execute in function for visualization of results
% %makeSegMovie(myelin, uint8(raw), '/home/alik/Desktop/test.avi');
% 
