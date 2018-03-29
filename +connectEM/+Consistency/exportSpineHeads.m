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
% NOTE(amotta): I've messed up `5115760af3200096c7847f016ce2dd3dd5b6a140`
% by not passing a seed value to `rng`. As a result, I cannot reproduce the
% spine heads used for the first set. Let's manually remove these spine
% heads based on their IDs (recovered from the NML file).
openShIds = [ ...
     70266, 304247, 269329, 187491, 227442, 123208, 309636,  78566, ...
    285553, 153211, 260123, 324407,  33731, 386422, 202399, 181221, ...
    185763, 127374, 211426, 212367, 330470];
openShIds = setdiff(1:numel(shAgglos), openShIds);

rng(0);
randIds = randperm(numel(openShIds));
randIds = openShIds(randIds(1:25));

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
