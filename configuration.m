% In this file parmeters for a pipeline run are set
% Note (!!!) that a second run with the same configuration will overwrite
% previous one
% Please set all parameters needed for your dataset

% Choose were to store result of the calculations, make sure you have WRITE access;
p.saveFolder = '/gaba/u/kostalr/data/documents/segTry/';
% Define ROI, make sure no black (surround) pixel are within this bbox
% This can be copied directly from webKNOSSOS bounding box field
p.bbox_wK = [1810, 2058, 398, 2916, 3399, 1150]; 
% Define directory and file prefix and voxel size for KNOSSOS hierachy with raw data
% READABLE to you on gaba
p.raw.root = '/gaba/u/kostalr/data/documents/ratL4stackKnoss/1/';
p.raw.prefix = '2016-03-12_ex4_RN_CC_S1BF_L4_mag1';
p.raw.voxelSize = [12 12 30];
% This segmentation parameter controls over- vs. undersegmentation,
% decrease for more smaller segments and vice versa
p.seg.threshold = .25;

% Do not change this, will add other parameters that usually need no
% modification
p = setParameterSettings(p);
