% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);

conn = connectEM.Connectome.load(param, connFile);

% Loading spine head agglomerates
shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

%% Count number of spine heads overlapping with axons
axonLUT = Agglo.buildLUT(maxSegId, conn.axons);

spineHeadsOverlapping = cellfun( ...
    @(segIds) any(axonLUT(segIds)), shAgglos);
numSpineHeadsOverlapping = sum(spineHeadsOverlapping) %#ok

spineHeadsContained = cellfun( ...
    @(segIds) all(axonLUT(segIds)), shAgglos);
numSpineHeadsContained = sum(spineHeadsContained) %#ok

%% Export contained spine heads
rng(0);
randIds = find(spineHeadsContained);
randIds = randIds(randperm(numel(randIds)));
randIds = randIds(1:25);

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

skel = Skeleton.fromMST(cellfun( ...
    @(segIds) segPoints(segIds, :), ...
    shAgglos(randIds), 'UniformOutput', false), ...
    param.raw.voxelSize, skel);
