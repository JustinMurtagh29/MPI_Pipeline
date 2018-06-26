% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

nmlDirs = '/mnt/mpibr/data/Personal/mottaa/L4/2018-05-13-Whole-Cell-Ground-Truth';
nmlDirs = fullfile(nmlDirs, { ...
    'HW_L4_WholeCellsCenterGTROI2017_13_10_2017'; ...
    'HW_L4_WholeCellsBorderGTROI2017_06_02_2018'});

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

%% Finding all NML files
nmlFiles = cellfun( ...
    @(nmlDir) dir(fullfile(nmlDir, '*.nml')), ...
    nmlDirs, 'UniformOutput', false);
nmlFiles = cat(1, nmlFiles{:});
nmlFiles = arrayfun( ...
    @(f) fullfile(f.folder, f.name), ...
    nmlFiles, 'UniformOutput', false);

%% Calculating path lengths
pathLensNm = cellfun(@(nmlFile) forNmlFile(param, nmlFile), nmlFiles);
pathLensTotalM = sum(pathLensNm) / 1E9;

%% Utilities
function pathLenNm = forNmlFile(param, nmlFile)
    pathLenNm = skeleton(nmlFile);
    pathLenNm = pathLenNm.pathLength([], param.raw.voxelSize);
    pathLenNm = sum(pathLenNm);
end
