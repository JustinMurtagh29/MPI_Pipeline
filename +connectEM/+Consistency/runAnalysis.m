% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

thisDir = fileparts(mfilename('fullpath'));
ctrlDir = fullfile(thisDir, 'annotations', 'random-spine-synapses');

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

%% Loading control data
ctrlFiles = dir(fullfile(ctrlDir, '*.nml'));
ctrlFiles = fullfile(ctrlDir, {ctrlFiles.name});

ctrlSynT = connectEM.Consistency.loadAnnotations(param, ctrlFiles);
ctrlSynT = vertcat(ctrlSynT{:});