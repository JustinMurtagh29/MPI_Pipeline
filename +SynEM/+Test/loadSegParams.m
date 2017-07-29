function p = loadSegParams( dataFolder )
%LOADSEGPARAMS Load the segmentation parameter struct for the SynEM test
%set.
% INPUT dataFolder: string
%           Path to SynEM data folder.
% OUTPUT p: struct
%           Segmentation parameter struct for test set.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

segFolder = fullfile(dataFolder, 'ex145_07x2', '20141007T094904');
rawFolder = fullfile(dataFolder, 'ex145_07x2', 'color', '1');
m = load(fullfile(segFolder, 'allParameter.mat'));
p = m.p;

%set path to test set cube
p = Seg.Util.changePaths(p, segFolder);

%paths to raw data and segmentation
p.raw.root =  rawFolder;
p.seg.root = fullfile(segFolder, 'globalSeg');

end

