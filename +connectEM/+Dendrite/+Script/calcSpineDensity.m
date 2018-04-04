% This script attempts to calculate the dendrite path length along the
% trunk. This is done by calculating the MST after removing spine head and
% neck segments.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
sansSpineFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_v2.mat');
withSpineFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_spine_attachment.mat');

info = Util.runInfo();

%% Loading data
Util.log('Loading data');
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);

sansSpine = load(sansSpineFile);
sansSpine = sansSpine.dendrites(sansSpine.indBigDends);
sansSpine = Agglo.fromSuperAgglo(sansSpine);

withSpine = load(withSpineFile);
shAgglos = withSpine.shAgglos;
withSpine = withSpine.dendAgglos(withSpine.indBigDends);

%% Remove spine heads and necks from dendrites
Util.log('Removing spine heads and necks');
mask = false(maxSegId, 1);
mask(cell2mat(sansSpine)) = true;
mask(cell2mat(shAgglos)) = false;

trunkAgglos = cellfun( ...
    @(segIds) segIds(mask(segIds)), ...
    withSpine, 'UniformOutput', false);

trunkSuperAgglos = cellfun( ...
    @(segIds) segPoints(segIds, :), ...
    trunkAgglos, 'UniformOutput', false);
trunkSuperAgglos = struct('nodes', trunkSuperAgglos);

%% Calculate path lengths
Util.log('Calculating path length');
trunkLensUm = Superagglos.mstLength( ...
    trunkSuperAgglos, param.raw.voxelSize);
trunkLensUm = trunkLensUm / 1E3;

%% Calculate spine density
Util.log('Calculating spine density');
shLUT = Agglo.buildLUT(maxSegId, shAgglos);
spineIds = cellfun( ...
    @(segIds) setdiff(shLUT(segIds), 0), ...
    withSpine, 'UniformOutput', false);

spineCount = cellfun(@numel, spineIds);
spineDensity = spineCount ./ trunkLensUm;

%% Plot spine densities
minPathLenUm = 10;
minPathMask = (trunkLensUm >= minPathLenUm);

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [480, 440];

ax = axes(fig);
axis(ax, 'square');
hold(ax, 'on');

binEdges = linspace(0, 1.5, 51);
histogram(ax, ...
    spineDensity, binEdges, ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);
histogram(ax, ...
    spineDensity(minPathMask), binEdges, ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);

ax.TickDir = 'out';
xlabel(ax, 'Spine density (µm^{-1})');
ylabel(ax, 'Dendrites');

legend(ax, ...
    'All dendrites', ...
    sprintf('All dendrites >= %d µm', minPathLenUm));

annotation(fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'HorizontalAlignment', 'center', ...
    'EdgeColor', 'none');
