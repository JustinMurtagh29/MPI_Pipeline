% Hi there!
%   Please do not remove the following lines of code.
%   If you cannot resist the temptation then godspeed!
if ~exist('PIPELINE_READY', 'var') || ~PIPELINE_READY
    error('Please start MATLAB inside the pipeline directory');
end
if ~strcmp(version('-release'), '2015b')
    error('Please run the pipeline code with Matlab R2015b, see README.md');
end

%%%%%%%%%%%%%%%%%%%%%%%%%% EDIT BELOW THIS LINE %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose were to store result of the calculations
% Make sure you have WRITE access
p.saveFolder = '/gaba/u/mberning/results/pipeline/20170217_ROI/';

% Define region of interest
% This can be copied directly from webKNOSSOS bounding box field
% Make sure p.bbox is always 25 pixels away from any black region in X-Y (50 in Z)
p.bbox_wK = [128, 128, 128, 5446, 8381, 3286]; 

% Name of the experiment. It's the same as the Dataset name on webKnossos.
% Also in the "Info" section when you open your dataset in webKnossos.
p.experimentName = '2012-09-28_ex145_07x2_ROI2017';

% Define directory and file prefix and voxel size for KNOSSOS hierachy
% with raw data READABLE to you on gaba
p.raw.root = '/gaba/u/mberning/wkCubes/2012-09-28_ex145_07x2_ROI2016_corrected/color/1/';
p.raw.prefix = '2012-09-28_ex145_07x2_ROI2016_corrected_mag1';

% Voxel size in nano metres
p.raw.voxelSize = [11.24 11.24 28];

% This segmentation parameter controls over- vs. undersegmentation,
% decrease for more smaller segments and vice versa
p.seg.threshold = .25;

% If p.myelin.isUsed is set to true a previously run myelin detection 
% (see preprocessing/additionalHeuristics.m) will be used to ensure that segments
% do not cross myelin/non-myelin border 
p.myelin.isUsed = true;
p.myelin.root = '/gaba/u/mberning/wkCubes/2012-09-28_ex145_07x2_ROI2016_vessel/segmentation/1/';
p.myelin.prefix = '2012-09-28_ex145_07x2_ROI2016_vessel_mag1';
p.myelin.segId = 3; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% STOP EDITING HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do not change this, will add other parameters that
% usually need no modification
p = setParameterSettings(p);

