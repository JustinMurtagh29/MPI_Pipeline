% This script assesses the errors introduced by the removal of overlaps
% between the axon and dendrite agglomerates.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connOldFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');
connNewFile = fullfile(rootDir, 'connectomeState', 'connectome_ax18a_deWC01wSp.mat');
eqClassFile = fullfile(rootDir, 'SVGDB', 'agglos', 'ax18a_deWC01wSp', 'eClass.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_spine_attachment.mat');
outDir = '/home/amotta/Desktop';

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);

connOld = load(connOldFile);
connNew = load(connNewFile);

eqClassLUT = load(eqClassFile);
eqClassLUT = eqClassLUT.segMapping;

shAgglos = load(shFile);
shAgglos = shAgglos.shAgglos;

%% Generate NML file with somata
somaAgglos = (connNew.denMeta.targetClass == 'Somata');
somaAgglos = connNew.dendrites(somaAgglos);

% Only look at ten somata
somaAgglos = somaAgglos(1:10);
somaCompIds = cell2mat(somaAgglos);

% Meta
skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

% Trees
skelNodes = arrayfun( ...
    @(id) segPoints(eqClassLUT == id, :), ...
    somaCompIds, 'UniformOutput', false);
skel = Skeleton.fromMST(skelNodes, param.raw.voxelSize, skel);

skel.write(fullfile(outDir, 'somata.nml'));

%% Check how many spine agglomerates vanished
% Find spine heads that are part of "old" dendrites
oldDendLUT = Agglo.buildLUT(maxSegId, connOld.dendrites);

newDendIds = cell2mat(connNew.dendrites);
[~, newDendLUT] = ismember(eqClassLUT, newDendIds);

oldMask = cellfun(@(ids) any(oldDendLUT(ids)), shAgglos);
newMask = cellfun(@(ids) any(newDendLUT(ids)), shAgglos);

vanishedMask = oldMask & (~newMask);
vanishedCount = sum(vanishedMask)