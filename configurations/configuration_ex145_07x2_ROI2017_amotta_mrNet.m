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
p.saveFolder = '/tmpscratch/amotta/l4/2018-10-10-mrnet-pipeline-run/';

% Define region of interest
% This can be copied directly from webKNOSSOS bounding box field
p.bbox_wK = [128, 128, 128, 5446, 8381, 3286];

% Name of the experiment. It's the same as the Dataset name on webKnossos.
% Also in the "Info" section when you open your dataset in webKnossos.
p.experimentName = '2012-09-28_ex145_07x2_ROI2017_MrNet_LogDistSeg';

% Define directory and file prefix and voxel size for KNOSSOS hierachy
% with raw data READABLE to you on gaba
p.raw = struct;
p.raw.root = '/gaba/u/mberning/wkCubes/2012-09-28_ex145_07x2_ROI2016_corrected/color/1';
p.raw.prefix = '2012-09-28_ex145_07x2_ROI2016_corrected_mag1';
p.raw.voxelSize = [11.24, 11.24, 28];
p.raw.dtype = 'uint8';

p.norm = struct;
p.norm.meanVal = 122.4522;
p.norm.stdVal = 20.5668;

p.class = struct;
p.class.root = '/tmpscratch/bstaffle/data/2012-09-28_ex145_07x2_ROI2017/mr_net_2e_cont5';
p.class.backend = 'wkwrap';
p.class.dtype = 'single';

% This segmentation parameter controls over- vs. undersegmentation,
% decrease for more smaller segments and vice versa
p.seg = struct;
p.seg.func = @(in) watershedSegMrNet( ...
    in, 'classThresh', 0, 'minDepth', 0.15, 'voxelSize', p.raw.voxelSize);
p.seg.backend = 'wkwrap';
p.seg.dtype = 'uint32';

%%%%%%%%%%%%%%%%%%%%%%%%%%%% STOP EDITING HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do not change this, will add other parameters that
% usually need no modification
p = setParameterSettings(p);
