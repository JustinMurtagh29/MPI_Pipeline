% This script assesses the errors introduced by the removal of overlaps
% between the axon and dendrite agglomerates.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_ax18a_deWC01wSp.mat');
outDir = '/home/amotta/Desktop';

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);

conn = load(connFile);

%% Generate NML file with somata
somaAgglos = (conn.denMeta.targetClass == 'Somata');
somaAgglos = conn.dendrites(somaAgglos);

% Meta
skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

% Trees
skelNodes = cellfun( ...
    @(ids) segPoints(ids, :), ...
    somaAgglos, 'UniformOutput', false);
skel = Skeleton.fromMST(skelNodes, param.raw.voxelSize, skel);

skel.write(fullfile(outDir, 'somata.nml'));