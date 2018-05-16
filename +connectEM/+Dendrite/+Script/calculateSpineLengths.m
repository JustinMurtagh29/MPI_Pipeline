% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
trunkFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3.mat');
dendFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto-and-manual.mat');

% For the origin of this magic constant, see
% https://gitlab.mpcdf.mpg.de/connectomics/amotta/blob/534c026acd534957b84e395d697ac48b3cc6a7ad/matlab/+L4/+Spine/buildDendriteTrunkMask.m
trunkMinSize = 10 ^ 5.5;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

trunks = load(trunkFile);
trunks = Agglo.fromSuperAgglo(trunks.dendrites(trunks.indBigDends));
trunks = trunks(Agglo.calculateVolume(param, trunks) > trunkMinSize);

dendrites = load(dendFile);
spineHeads = dendrites.shAgglos;
dendrites = dendrites.dendrites(dendrites.indBigDends);

%% Calculate spine lengths
spineLengths = ...
    connectEM.Dendrite.calculateSpineLengths( ...
        param, trunks, dendrites, spineHeads);
