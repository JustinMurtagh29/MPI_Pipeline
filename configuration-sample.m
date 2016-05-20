% Hi there!
%   Please do not remove the follwing lines of code.
%   If you cannot resist the temptation then godspeed!
if ~exist('PIPELINE_READY', 'var') || ~PIPELINE_READY
    error('Please start MATLAB inside the pipeline directory');
end

%%%%%%%%%%%%%%%%%%%%%%%%%% EDIT BELOW THIS LINE %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose were to store result of the calculations
% Make sure you have WRITE access
p.saveFolder = '/gaba/u/mberning/results/pipeline/test/';

% Define ROI, make sure no black (surround) pixel are within this bbox
% This can be copied directly from webKNOSSOS bounding box field
p.bbox_wK = [1153, 769, 129, 7808, 5376, 3200]; 
  
% Define directory and file prefix and voxel size for KNOSSOS hierachy
% with raw data READABLE to you on gaba
p.raw.root = '/gaba/u/mberning/data/cortex/2012-09-28_ex145_07x2_corrected/color/1/';
p.raw.prefix = '2012-09-28_ex145_07x2_corrected_mag1';

% Voxel size in nano metres
p.raw.voxelSize = [11.24 11.24 28];

% This segmentation parameter controls over- vs. undersegmentation,
% decrease for more smaller segments and vice versa
p.seg.threshold = .25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% STOP EDITING HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do not change this, will add other parameters that
% usually need no modification
p = setParameterSettings(p);