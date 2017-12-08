% Write out location of random spine heads, which then serve as seeds for
% groundtruth dendrite (center-line) tracings.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

%% load data
paramFile = fullfile(rootDir, 'allParameter.mat');
param = load(paramFile, 'p');
param = param.p;

postSynFile = fullfile( ...
    rootDir, 'aggloState', ...
    'dend+wholecells_01_auto.mat');
load(postSynFile, 'shAgglos');

points = Seg.Global.getSegToPointMap(param);

%% select random spine heads
numSpineHeads = 50;

rng(0);
randIds = randperm(numel(shAgglos));
randIds = randIds(1:numSpineHeads);

randShAgglos = shAgglos(randIds(:));

%% build skeleton
randShPoints = cell2mat(cellfun( ...
    @(segIds) points(segIds(1), :), ...
    randShAgglos, 'UniformOutput', false));

skel = skeleton();
skel = skel.addNodesAsTrees(randShPoints);

skel = Skeleton.setParams4Pipeline(skel, param);
skel.write('/home/amotta/Desktop/random-spine-heads.nml');