% This script helps to analyse some of the differences found between the
% connectome versions `connectome_axons_18_a_ax_spine_syn_clust` and
% `connectome_ax18a_deWC01wSp` (i.e., the edel connectome).
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

connOldFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');
synOldFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');

connNewFile = fullfile(rootDir, 'connectomeState', 'connectome_ax18a_deWC01wSp.mat');
synNewFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v4_ax18a_deWC01wSp.mat');

edgesNewFile = fullfile(rootDir, 'SVGDB', 'agglos', 'ax18a_deWC01wSp', 'edges.mat');
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
synOld = load(synOldFile);

connNew = load(connNewFile);
synNew = load(synNewFile);

eqClassLUT = load(eqClassFile);
eqClassLUT = eqClassLUT.segMapping;

edgesNew = load(edgesNewFile);

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
vanishedCount = sum(vanishedMask);

%% Check how many spine synapses vanished
% Compute expected spine synapses in "new" segmentation.
%
% NOTE(amotta): For some reason the indexing into `borderMappin` produces
% an index error. The maximum `edgeIdx` is larger than the rows in:
% * `edges` of `/gaba/u/mberning/results/pipeline/20170217_ROI/globalEdges.mat`
% * `borderMapping` of `/gaba/u/mberning/results/pipeline/20170217_ROI/SVGDB/agglos/ax18a_deWC01wSp/edges.mat`
newSpineSynEdgesExp = ...
    synOld.synapses.edgeIdx(synOld.isSpineSyn);
newSpineSynEdgesExp = cellfun( ...
    @(ids) edgesNew.borderMapping(ids), ...
    newSpineSynEdgesExp, 'UniformOutput', false);

% Compare against actually found spine synapses
newSpineSynEdges = synNew.synapses.edgeIdx(synNew.isSpineSyn);
newSpineSynEdges = cell2mat(newSpineSynEdges);

%% Plot distribution of spine fraction
% This is important for the definition of axon classes.
connOld.axonMeta.spineSynFrac = ...
    connOld.axonMeta.spineSynCount ...
 ./ connOld.axonMeta.synCount;
maskOld = (connOld.axonMeta.synCount >= 10);

connNew.axonMeta.spineSynFrac = ...
    connNew.axonMeta.spineSynCount ...
 ./ connNew.axonMeta.synCount;
maskNew = (connNew.axonMeta.synCount >= 10);

% Plotting
binEdges = linspace(0, 1, 21);

fig = figure();
fig.Color = 'white';

ax = axes(fig);
axis(ax, 'square');
hold(ax, 'on');

histogram( ...
    ax, connOld.axonMeta.spineSynFrac(maskOld), ...
    'BinEdges', binEdges, ...
    'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);
histogram( ...
    ax, connNew.axonMeta.spineSynFrac(maskNew), ...
    'BinEdges', binEdges, ...
    'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);

xlim(ax, binEdges([1, end]));
legend(ax, 'Old', 'Edel', 'Location', 'NorthWest');
xlabel('Fraction of synapses onto spines');
ylabel('Fraction of axons');

title( ...
    ax, {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Compare spine synapses sizes
oldSynT = table;
oldSynT.id = cell2mat(connOld.connectome.synIdx);
oldSynT.area = cell2mat(connOld.connectomeMeta.contactArea);
oldSynT.isSpine = synOld.isSpineSyn(oldSynT.id);

newSynT = table;
newSynT.id = cell2mat(connNew.connectome.synIdx);
newSynT.area = cell2mat(connNew.connectomeMeta.contactArea);
newSynT.isSpine = synNew.isSpineSyn(newSynT.id);

% Restrict to spine synapses
oldSynT(~oldSynT.isSpine, :) = [];
[~, uniRows] = unique(oldSynT.id);
oldSynT = oldSynT(uniRows, :);

newSynT(~newSynT.isSpine, :) = [];
[~, uniRows] = unique(newSynT.id);
newSynT = newSynT(uniRows, :);

% Plot
binEdges = [oldSynT.area; newSynT.area];
binEdges = linspace(0, prctile(binEdges, 95), 51);

fig = figure();
fig.Color = 'white';

ax = axes(fig);
ax.TickDir = 'out';
axis(ax, 'square');
hold(ax, 'on');

plotIt = @(data) ...
    histogram( ...
        ax, data, ...
        'BinEdges', binEdges, ...
        'Normalization', 'probability', ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
plotIt(oldSynT.area);
plotIt(newSynT.area);

xlim(ax, binEdges([1, end]));
xlabel(ax, 'Spine synapse area (µm²)');
ylabel(ax, 'Fraction of spine synapses');

legend(ax, ...
    'SegEM-based connectome', ...
    'Edel connectome', ...
    'Location', 'NorthEast');
title(ax, ...
   {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
