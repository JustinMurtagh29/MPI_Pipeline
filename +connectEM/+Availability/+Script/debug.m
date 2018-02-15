% The aim of this script is to check whether the results of the
% availability calculations do make sense and to find remaining issues.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
outDir = '/home/amotta/Desktop';
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
availFile = '/tmpscratch/amotta/l4/2018-02-02-surface-availability-connectome-axons-18-a/axon-avail-data.mat';

connFile = fullfile( ...
    rootDir, 'connectomeState', ...
    'connectome_axons_18_a_with_den_meta.mat');

className = 'ApicalDendrite';
radiusNm = 1000;

minSpecificity = 0.2;
minSynCount = 10;

info = Util.runInfo();

%% loading data
avail = load(availFile);
conn = load(connFile);

param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

points = Seg.Global.getSegToPointMap(param);

%%
avail.axonAvail = avail.axonAvail(2:end, :, :);
avail.targetClasses = avail.targetClasses(2:end);

%% find axons with smallest & largest availabilities
axonIds = find(conn.axonMeta.synCount >= minSynCount);

axonAvail = squeeze(avail.axonAvail( ...
    :, avail.dists == radiusNm, axonIds));
axonAvail = axonAvail ./ sum(axonAvail, 1);
axonAvail = axonAvail( ...
    avail.targetClasses == className, :);

% sort by increasing availability
axons = table;
[axons.avail, axons.id] = sort(axonAvail(:), 'ascend');
axons.id = axonIds(axons.id);

%% generate skeleton
skel = skeleton();

postAgglos = (conn.denMeta.targetClass == className);
postAgglos = conn.dendrites(postAgglos);

skel = Skeleton.fromMST( ...
    cellfun( ...
        @(segIds) points(segIds, :), ...
        postAgglos, 'UniformOutput', false), ...
	param.raw.voxelSize, skel);
skel.names(:) = {className};
skel.colors(:) = {[0, 0, 1, 1]};

preAxons = axons([1:10, (end - 9):end], :);
preAgglos = conn.axons(preAxons.id);

preNames = arrayfun( ...
    @(id, avail) sprintf('Axon %d, Availability %.2f', id, avail), ...
    preAxons.id, preAxons.avail, 'UniformOutput', false);

skel = Skeleton.fromMST( ...
    cellfun( ...
        @(segIds) points(segIds, :), ...
        preAgglos, 'UniformOutput', false), ...
	param.raw.voxelSize, skel);
skel.names((end - numel(preAgglos) + 1):end) = preNames;
skel.colors((end - numel(preAgglos)     + 1):end) = {[1, 0, 0, 1]};
skel.colors((end - numel(preAgglos) / 2 + 1):end) = {[0, 1, 0, 1]};

skel = Skeleton.setParams4Pipeline(skel, param);

skelFile = sprintf('%s-avail-debug.nml', lower(className));
skel.write(fullfile(outDir, skelFile));