% Hi there!
%   Please do not remove the following lines of code.
%   If you cannot resist the temptation then godspeed!
global PIPELINE_READY
if isempty(PIPELINE_READY)
    error('Please start MATLAB inside the pipeline directory');
end
clear PIPELINE_READY

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
p.saveFolder = '/tmpscratch/kunsong/EE_M2_6R_SegEMTest2_tiltCorrected_seg2_20181019/';

% Define region of interest
% This can be copied directly from webKNOSSOS bounding box field
% Make sure p.bbox is always 25 pixels away from any black region in X-Y (10 in Z)
p.bbox_wK = [256,256,128,5800,7610,5362];

% Name of the experiment. It's the same as the Dataset name on webKnossos.
% Also in the "Info" section when you open your dataset in webKnossos.
p.experimentName = 'KS_EE_round1_20160810_tiltCorrected';

% If you want to be notified via email after the completion of your pipeline run
p.email.notify = false; %default no notification
if p.email.notify
    p.email.address = 'kun.song@brain.mpg.de';
end
  
% Define directory and file prefix and voxel size for KNOSSOS hierachy
% with raw data READABLE to you on gaba
p.raw.root = '/tmpscratch/kunsong/EE_M2_6R_SegEMTest2_tiltCorrected_wkw/color/1/';
% p.raw.prefix = 'KS_EE_round1_20160810_tiltCorrected_mag1';

% Uncomment and change this, if your dataset comprises a mask knossos hierarchy (e.g. created
% by the KSMB cubing), which marks the outside (e.g. padded) region of the
% dataset
% p.mask.root = '/gaba/u/mberning/data/cortex/2012-09-28_ex145_07x2_corrected/mask/1/'
% p.mask.prefix = '2012-09-28_ex145_07x2_corrected_mag1_mask';

% Uncomment, if you wanna live the risky life
p.raw.backend = 'wkwrap';

% Voxel size in nano metres
p.raw.voxelSize = [12 12 30];

% This segmentation parameter controls over- vs. undersegmentation,
% decrease for more smaller segments and vice versa
p.seg.threshold = .23;

% If p.myelin.isUsed is set to true a previously run myelin detection 
% (see preprocessing/additionalHeuristics.m) will be used to ensure that segments
% do not cross myelin/non-myelin border 
p.myelin.isUsed = false;
p.myelin.root = '/gaba/u/mberning/wkCubes/2012-09-28_ex145_07x2_ROI2016_vessel/segmentation/1/';
p.myelin.prefix = '2012-09-28_ex145_07x2_ROI2016_vessel_mag1';
p.myelin.segId = 3; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% STOP EDITING HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do not change this, will add other parameters that
% usually need no modification
p = setParameterSettings(p);

