% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

wkParam = struct;
wkParam.root = '/gaba/wKcubes/Connectomics department/2012-09-28_ex145_07x2_ROI2017/color/1/';
wkParam.prefix = '2012-09-28_ex145_07x2_ROI2016_mag1';

roiSize = 1.5 .* [512; 512; 256] + 1;
padSize = [25; 25; 10];

wkOffset = [129; 129; 129];
wkwOffset = [513; 513; 513];

%% build base configuration
rootDir = fullfile('/tmpscratch', getenv('USER'));
rootDir = fullfile(rootDir, 'wkwrap', datestr(now, 30));
assert(not(exist(rootDir, 'dir')));
mkdir(rootDir);

% Choose were to store result of the calculations
% Make sure you have WRITE access
p = struct;
p.saveFolder = fullfile(rootDir, 'pipeline');
mkdir(p.saveFolder);

% Define region of interest
% This can be copied directly from webKNOSSOS bounding box field
% Make sure p.bbox is always 25 pixels away from any black region in X-Y (10 in Z)
p.bbox_wK = [wkwOffset - 1, roiSize];
p.bbox_wK = reshape(p.bbox_wK, 1, []);

% Name of the experiment. It's the same as the Dataset name on webKnossos.
% Also in the "Info" section when you open your dataset in webKnossos.
p.experimentName = 'wkwrap-test';
  
% Define directory and file prefix and voxel size for KNOSSOS hierachy
% with raw data READABLE to you on gaba
p.raw.root = fullfile(rootDir, 'wkw', 'color', '1');
p.raw.backend = 'wkwrap';

% Voxel size in nano metres
p.raw.voxelSize = [11.24 11.24 28];

% This segmentation parameter controls over- vs. undersegmentation,
% decrease for more smaller segments and vice versa
p.seg.threshold = .25;

% If p.myelin.isUsed is set to true a previously run myelin detection 
% (see preprocessing/additionalHeuristics.m) will be used to ensure that segments
% do not cross myelin/non-myelin border 
p.myelin.isUsed = false;

%% copy ROI from KNOSSOS to .wkw
copyBoxWk = bsxfun(@plus, wkOffset, [zeros(3, 1), roiSize - 1]);
copyBoxWk = copyBoxWk + bsxfun(@times, padSize, [-1, +1]);
copyOffWkw = wkwOffset - padSize;

fprintf('Copying raw data... ');
raw = readKnossosRoi(wkParam.root, wkParam.prefix, copyBoxWk);
wkwInit('new', p.raw.root, 32, 32, 'uint8', 1);
saveRawData(p.raw, copyOffWkw, raw);
disp('Done!');

%% complete parameter structure
p = setParameterSettings(p);

%% run pipeline (and see how far we get...)
runPipeline(p);