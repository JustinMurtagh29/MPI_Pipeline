% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
somaFile = fullfile(rootDir, 'aggloState', 'somata_06.mat');
dendFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v2_auto-and-manual.mat');

info = Util.runInfo();

%% Loading data
Util.log('Loading data');

param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

soma = load(somaFile);
dend = load(dendFile);

%% Remove somata from dendrite super-agglomerates
Util.log('Removing somata from dendrites');

somaSegIds = cell2mat(Agglo.fromSuperAgglo(soma.somata));
dendrites = SuperAgglo.removeSegIds(param, dend.dendrites, somaSegIds);

% TODO(amotta): Enable check if possible
dendrites = SuperAgglo.clean(dendrites, false);

%% Add somata to list of postsynaptic processes
% TODO(amotta): Compute whole-cell indices

%% Completing and saving result
% TODO(amotta): Remove empty super-agglomerates
