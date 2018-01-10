% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_16_b.mat');
outDir = '/tmpscratch/amotta/l4/2018-01-10-all-axons-for-amira';

%% loading data
fprintf('Loading data... ');
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

axons = load(axonFile, 'axons', 'indBigAxons');
axons = axons.axons(axons.indBigAxons);
fprintf('done!\n');

%% exporting to amira
Superagglos.exportToAmira(param, axons, outDir);