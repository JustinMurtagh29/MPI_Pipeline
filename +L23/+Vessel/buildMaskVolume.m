% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/tmpscratch/amotta/l23/2018-10-09-mrnet-pipeline-run';
outDir = '/tmpscratch/amotta/l23/2018-10-13-vessel-and-nuclei-mask';

% Segmentation for further annotation in webKnossos
annSegOut = struct;
annSegOut.root = '/tmpscratch/webknossos/Connectomics_Department/2018-10-13_ex144_08x2_blood-vessel-mask/segmentation/1';
annSegOut.backend = 'wkwrap';

mag = [8, 8, 4];

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading parameters
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

%% Loading segmentation
rawParam = param.raw;
rawParam.root = fullfile( ...
    fileparts(rawParam.root(1:(end - 1))), ...
    sprintf('%d-%d-%d', mag));

magBox = ceil(param.bbox ./ mag(:));
raw = loadRawData(rawParam, magBox);

%% Build mask
mask = raw > 198;
mask = bwareaopen(mask, 540000);
mask = imclose(mask, strel('cube', 5));

%% Export to webKnossos for further annotation
if exist('annSegOut', 'var')
    wkwInit('new', annSegOut.root, 32, 32, 'uint32', 1);
    wkwSaveRoi(annSegOut.root, magBox(:, 1)', uint32(mask));
    % NOTE(amotta): color/1 is symbolic link to `rawParam.root`.
end
