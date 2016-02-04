
% Choose were to store result of the calculations, make sure you have WRITE access;
p.saveFolder = '/gaba/u/mberning/results/pipeline/test/';
% Define ROI, numbers should be aligned with KNOSSOS cubes
p.bbox = [1153 7808; 769 5376; 129 3200]; 
% Define directory and file prefix for KNOSSOS hierachy with raw data
p.raw.root = '/gaba/u/mberning/data/cortex/2012-09-28_ex145_07x2_corrected/color/1/';
p.raw.prefix = '2012-09-28_ex145_07x2_corrected_mag1';

% Add constant and optional settings to the parameter structure
p = setParameterSettings(p);
% This will take some time, run classification, segmentation, correspondences, globalization, graph construction
runPipeline(p);

