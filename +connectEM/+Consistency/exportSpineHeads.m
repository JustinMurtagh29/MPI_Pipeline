% Exports random spine heads for manual axon-spine interface annotations.
% This data will be used for the control condition.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_spine_attachment.mat');
outputDir = '/home/amotta/Desktop';

% Number of spine heads to export
shCount = 25;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);

shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

%% Export random examples
rng(0);
randIds = randperm(numel(shAgglos), 25);
randShAgglos = cellfun( ...
    @(segIds) segPoints(segIds, :), ...
    shAgglos(randIds), 'UniformOutput', false);

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));
skel = Skeleton.fromMST(randShAgglos, param.raw.voxelSize, skel);

skel.names = arrayfun( ...
    @(idx, shId) sprintf('Syn %d, Spine head %d', idx, shId), ...
    1:numel(randIds), randIds, 'UniformOutput', false);
skel.write(fullfile(outputDir, 'rand-spine-heads.nml'));
