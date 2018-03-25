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

connNewFile = fullfile(rootDir, 'connectomeState', 'connectome_ax18a_deWC01wSp_v4.mat');
synNewFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v5_ax18a_deWC01wSp_withShId.mat');

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

%% Quantitative comparison of synapses
syns = {synOldFile, synNewFile};
data = cell(0, numel(syns) + 1);

for curIdx = 1:numel(syns)
    curSyn = syns{curIdx};
    curSyn = load(curSyn);
    
    data{1, 1} = '# synapses';
    data{1, 1 + curIdx} = numel(curSyn.isSpineSyn);
    data{2, 1} = '# spine synapses';
    data{2, 1 + curIdx} = sum(curSyn.isSpineSyn);
end

t = cell2table(data(:, 2:end));
t.Properties.RowNames = data(:, 1);
[~, t.Properties.VariableNames] = cellfun( ...
    @fileparts, syns, 'UniformOutput', false);
disp(t)

%% Look at number of spine heads with synapses
shLUT = Agglo.buildLUT(maxSegId, shAgglos);
synOld.ontoSpineHeadId = cellfun( ...
    @(segIds) max(shLUT(segIds)), ...
    synOld.synapses.postsynId);

%% Compare spine heads between connectomes
shIdsConnOld = cell2mat(connOld.connectome.synIdx);
shIdsConnOld = nonzeros(synOld.ontoSpineHeadId(shIdsConnOld));
fprintf('# spine heads with synapse (old): %d\n', numel(unique(shIdsConnOld)));

shIdsConnNew = cell2mat(connNew.connectome.synIdx);
shIdsConnNew = nonzeros(synNew.ontoSpineHeadId(shIdsConnNew));
fprintf('# spine heads with synapse (new): %d\n', numel(unique(shIdsConnNew)));

fprintf('\n');
fprintf('# spine heads lost: %d\n', numel(setdiff(shIdsConnOld, shIdsConnNew)));
fprintf('# spine heads added: %d\n', numel(setdiff(shIdsConnNew, shIdsConnOld)));

%% Check for split interfaces
connOld.axonMeta.synIds = accumarray( ...
    connOld.connectome.edges(:, 1), ...
    transpose(1:size(connOld.connectome, 1)), size(connOld.axons), ...
    @(ids) {cell2mat(connOld.connectome.synIdx(ids))}, {zeros(0, 1)});
connOld.axonMeta.shIds = cellfun( ...
    @(synIds) setdiff(synOld.ontoSpineHeadId(synIds), 0), ...
    connOld.axonMeta.synIds, 'UniformOutput', false);
splitShSynOld = sum( ...
    connOld.axonMeta.spineSynCount ...
  - cellfun(@numel, connOld.axonMeta.shIds));

connNew.axonMeta.synIds = accumarray( ...
    connNew.connectome.edges(:, 1), ...
    transpose(1:size(connNew.connectome, 1)), size(connNew.axons), ...
    @(ids) {cell2mat(connNew.connectome.synIdx(ids))}, {zeros(0, 1)});
connNew.axonMeta.shIds = cellfun( ...
    @(synIds) setdiff(synNew.ontoSpineHeadId(synIds), 0), ...
    connNew.axonMeta.synIds, 'UniformOutput', false);
splitShSynNew = sum( ...
    connNew.axonMeta.spineSynCount ...
  - cellfun(@numel, connNew.axonMeta.shIds));

fprintf('\n');
fprintf('# split ASIs (old): %d\n', splitShSynOld);
fprintf('# split ASIs (new): %d\n', splitShSynNew);

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
binEdges = linspace(0, prctile(binEdges, 99), 51);

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
