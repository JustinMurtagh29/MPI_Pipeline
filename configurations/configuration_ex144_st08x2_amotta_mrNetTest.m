% Hi there!
%   Please do not remove the following lines of code.
%   If you cannot resist the temptation then godspeed!
global PIPELINE_READY

if isempty(PIPELINE_READY)
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
p.saveFolder = '/tmpscratch/amotta/l23/2018-10-09-mrnet-segmentation-test/';

% Define region of interest
% This can be copied directly from webKNOSSOS bounding box field
p.bbox_wK = [512, 512, 512, 3 * 512, 3 * 512, 3 * 256];

% Name of the experiment. It's the same as the Dataset name on webKnossos.
% Also in the "Info" section when you open your dataset in webKnossos.
p.experimentName = '2018-10-01_ex144_st08x2_mrnet-test';

% Define directory and file prefix and voxel size for KNOSSOS hierachy
% with raw data READABLE to you on gaba
p.raw = struct;
p.raw.root = '/tmpscratch/amotta/l23/2018-10-04-gradient-corrected/color/1';
p.raw.backend = 'wkwrap';
p.raw.voxelSize = [11.24, 11.24, 28];

p.norm = struct;
p.norm.meanVal = 165.2;
p.norm.stdVal = 19.8;

p.class = struct;
p.class.root = '/tmpscratch/amotta/l23/2018-10-06-pipeline-run/mrnet';
p.class.backend = 'wkwrap';
p.class.dtype = 'single';

% This segmentation parameter controls over- vs. undersegmentation,
% decrease for more smaller segments and vice versa
p.seg = struct;
p.seg.func = @(in) watershedSegMrNet( ...
    in, 'classThresh', 0, 'minDepth', 0.2, 'voxelSize', p.raw.voxelSize);
p.seg.backend = 'wkwrap';
p.seg.dtype = 'uint32';

%%%%%%%%%%%%%%%%%%%%%%%%%%%% STOP EDITING HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do not change this, will add other parameters that
% usually need no modification
p = setParameterSettings(p);
