% Hi there!
%   Please do not remove the following lines of code.
%   If you cannot resist the temptation then godspeed!
if ~exist('PIPELINE_READY', 'var') || ~PIPELINE_READY
    error('Please start MATLAB inside the pipeline directory');
end
if ~strcmp(version('-release'), '2015b')
    error('Please run the pipeline code with Matlab R2015b, see README.md');
end
if exist('p', 'var')
    error(['Variable ''p'' already exists and might cause a corrupted ' ...
        'segmentation parameter file. Please delete the ''p'' variable' ...
        ' before running this script.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%% EDIT BELOW THIS LINE %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose were to store result of the calculations
% Make sure you have WRITE access
p.saveFolder = '/tmpscratch/amotta/pipeline/20170808-ex145-07x2-beneNet/';

% Define region of interest
% This can be copied directly from webKNOSSOS bounding box field
% Make sure p.bbox is always 25 pixels away from any black region in X-Y (10 in Z)
p.bbox_wK = [128, 128, 128, 5446, 8381, 3286]; 

% Name of the experiment. It's the same as the Dataset name on webKnossos.
% Also in the "Info" section when you open your dataset in webKnossos.
p.experimentName = '20170808-ex145-07x2-beneNet';
  
% Define directory and file prefix and voxel size for KNOSSOS hierachy
% with raw data READABLE to you on gaba
p.raw.root = '/gaba/wKcubes/Connectomics department/2012-09-28_ex145_07x2_ROI2016_corrected/color/1/';
p.raw.prefix = '2012-09-28_ex145_07x2_ROI2016_corrected_mag1';

% Uncomment and change this, if your dataset comprises a mask knossos hierarchy (e.g. created
% by the KSMB cubing), which marks the outside (e.g. padded) region of the
% dataset
% p.mask.root = '/gaba/u/mberning/data/cortex/2012-09-28_ex145_07x2_corrected/mask/1/'
% p.mask.prefix = '2012-09-28_ex145_07x2_corrected_mag1_mask';

% Uncomment, if you wanna live the risky life
% p.raw.backend = 'wkwrap';

% Voxel size in nano metres
p.raw.voxelSize = [11.24 11.24 28];

% This segmentation parameter controls over- vs. undersegmentation,
% decrease for more smaller segments and vice versa
p.seg.threshold = 0.25;

% If p.myelin.isUsed is set to true a previously run myelin detection 
% (see preprocessing/additionalHeuristics.m) will be used to ensure that segments
% do not cross myelin/non-myelin border 
p.myelin.isUsed = true;
p.myelin.root = '/gaba/wKcubes/Connectomics department/2012-09-28_ex145_07x2_ROI2016_vessel/segmentation/1/';
p.myelin.prefix = '2012-09-28_ex145_07x2_ROI2016_vessel_mag1';
p.myelin.segId = 3; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% STOP EDITING HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do not change this, will add other parameters that
% usually need no modification
p = setParameterSettings(p);

