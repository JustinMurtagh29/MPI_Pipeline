%% In this file three examples use cases of the pipeline are presented
% Each of those is delinated by a line that starts with %% (matlab cell
% mode)
% First Example: Run classification and segmentation on small subset of
% data for testing (writing a movie to look at results)
% Second Example: Run on whole dataset to upload to webKNOSSOS & isosurface
% visualizations based on skeletons
% Third Example: Run synapse detection on whole dataset

%% O-th section: Set all parameters needed for your dataset
% Note that some of those will be changed again in the specific examples
% below

% Choose were to store result of the calculations, make sure you have WRITE access;
p.saveFolder = '/gaba/u/mberning/results/pipeline/test/';
% Define ROI, make sure no black (surround) pixel are within this bbox
% This can be copied directly from webKNOSSOS bounding box field
p.bbox_wK = [1153, 769, 129, 7808, 5376, 3200]; 
% Define directory and file prefix for KNOSSOS hierachy with raw data
% READABLE to you on gaba
p.raw.root = '/gaba/u/mberning/data/cortex/2012-09-28_ex145_07x2_corrected/color/1/';
p.raw.prefix = '2012-09-28_ex145_07x2_corrected_mag1';
p.raw.voxelSize = [11.24 11.24 28];
% This segmentation parameter controls over- vs. undersegmentation,
% decrease for more smaller segments and vice versa
p.seg.threshold.hmin = .25;

%% 1st example: Generate segmentation movie for small subset of data

% Change bounding box to some small region you would want to see a movie of
% This will provide 800 by 600 by 128 frame movie
% If you choose this numbers too large you might get memory errors
p.bbox = [1001 1600; 1001 1800; 129 256]; 
p.seg.threshold.hmin = .25;
% Just make a short movie for judgement of segmentation quality
outputFile = makeSegmentationPreviewMovie(p);
% Display name of movie just generated
display(outputFile);

%% 2nd example: Example parameters pipeline isosurface visualization

% Add constant and optional settings to the parameter structure
p.bbox = [1153 7808; 769 5376; 129 3200]; 
p.seg.threshold.hmin = .25;
% This will take some time, run classification, segmentation, correspondences, globalization, graph construction, feature calculation
runPipeline(setParameterSettings(p));

%% 3rd example: Example parameters pipeline synapse detection

% Add constant and optional settings to the parameter structure
p.bbox = [1153 7808; 769 5376; 129 3200]; 
p.seg.threshold.hmin = .39;
% This will take some time, run classification, segmentation, correspondences, globalization, graph construction, feature calculation
runPipeline(setParameterSettings(p));
