% Hi there!
%   Please do not remove the following lines of code.
%   If you cannot resist the temptation then godspeed!
global PIPELINE_READY
if isempty(PIPELINE_READY)
    error('Please start MATLAB inside the pipeline directory');
end
clear PIPELINE_READY

% if ~strcmp(version('-release'), '2015b')
%     error('Please run the pipeline code with Matlab R2015b, see README.md');
% end
if exist('p', 'var')
    error(['Variable ''p'' already exists and might cause a corrupted ' ...
        'segmentation parameter file. Please delete the ''p'' variable' ...
        ' before running this script.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%% EDIT BELOW THIS LINE %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose were to store result of the calculations
% Make sure you have WRITE access
p.saveFolder = '/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet_HC/';

% Define region of interest
% This can be copied directly from webKNOSSOS bounding box field
% Make sure p.bbox is always 25 pixels away from any black region in X-Y (10 in Z)
p.bbox_wK = [3232, 3232, 32,11136, 9472, 3584]; 

% Name of the experiment. It's the same as the Dataset name on webKnossos.
% Also in the "Info" section when you open your dataset in webKnossos.
p.experimentName = 'H2_3_v2_U1_SubI';
% If you want to be notified via email after the completion of your pipeline run
p.email.notify = true; %default no notification
if p.email.notify
    p.email.address = 'sahil.loomba@brain.mpg.de';
end

% Define directory and file prefix and voxel size for KNOSSOS hierachy
% with raw data READABLE to you on gaba
p.raw.root = '/tmpscratch/jast/wKcubes/H2_3_v2_U1_SubI/color/1/';
% p.raw.prefix = '2012-09-28_ex145_07x2_corrected_mag1';
p.raw.backend = 'wkwrap';

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
p.seg.threshold = 0.03;

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


%%%%%%%%%%%%%%%%%%%%% OVERWRITE SOME DEFAULT SETTINGS %%%%%%%%%%%%%%%%%%%%%

% Use same normalization values as in H2_3_v2_U1_SubI run
p.norm = struct;
p.norm.meanVal = 127.7505;
p.norm.stdVal = 48.0068;
p.norm.func = @(x) normalizeStack(x, p.norm.meanVal, p.norm.stdVal);

% CNN was applied by Benedikt
p = rmfield(p, 'cnn');

% Use output from Benedikt's CNN
p.class = struct;
p.class.root = '/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/mr_2e_cont5_slurm/';
p.class.backend = 'wkwrap';
p.class.dtype = 'single';

p.seg = struct;
p.seg.func = @(in) watershedSegMrNet( ...
    in, 'classThresh', 0, 'minDepth', 0.15, 'voxelSize', p.raw.voxelSize);
p.seg.backend = 'wkwrap';
p.seg.dtype = 'uint32';
p.seg.root = fullfile(p.saveFolder, 'globalSeg', '1/');
 

outFile = fullfile(p.saveFolder, 'allParameter.mat');
Util.save(outFile, p);


%%%%%%%%%%%%%%%%% START PIPELINE (WITHOUT CLASSIFICATION) %%%%%%%%%%%%%%%%%

runPipeline(p, PipelineStep.MyelinFix, PipelineStep.CompressSegmentation);
%{
% 21.09.2018
% for documentation (this was run manually and not from this script)
% not sure why I skipped that in the first place
runPipeline(p, PipelineStep.GraphConstruction, ...
    PipelineStep.GlobalGraphStruct);
%}
