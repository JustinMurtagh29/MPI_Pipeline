% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
autoFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');
fullFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto-and-manual.mat');

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

auto = load(autoFile);
full = load(fullFile);

%% Show results
datasetVol = 1 + param.bbox(:, 2) - param.bbox(:, 1);
datasetVol = prod(param.raw.voxelSize(:) / 1E3 .* datasetVol);

numSpineHeads = numel(auto.shAgglos)
spineHeadDensity = numSpineHeads / datasetVol
numAutoAttached = sum(auto.attached ~= 0)
fracAutoAttached = numAutoAttached / numSpineHeads
numManualAttached = sum(full.attached & ~auto.attached)

% TODO(amotta): Calculate automatically obtained spine density
