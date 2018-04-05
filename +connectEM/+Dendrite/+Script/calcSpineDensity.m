% This script attempts to calculate the dendrite path length along the
% trunk. This is done by calculating the MST after removing spine head and
% neck segments.
%
% TODO
% * Exclude whole cells by brute force, or
% * Make sure that only IN whole cells are detected
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
trunkFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_v2.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_spine_attachment.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');
outDir = '/home/amotta/Desktop';

info = Util.runInfo();

%% Loading data
Util.log('Loading data');
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);

% Dendrite prior to spine attachment
trunks = load(trunkFile);
trunks = trunks.dendrites(trunks.indBigDends);
trunks = Agglo.fromSuperAgglo(trunks);

% Spine heads
shAgglos = load(shFile);
shAgglos = shAgglos.shAgglos;

% Connectome (and dendrites after spine attachment)
conn = load(connFile);
dendrites = conn.dendrites;

somata = (conn.denMeta.targetClass == 'Somata');
somata = conn.dendrites(somata);

%% Remove spine heads and necks from dendrites
Util.log('Removing spine heads and necks');
mask = false(maxSegId, 1);
mask(cell2mat(trunks)) = true;
mask(cell2mat(shAgglos)) = false;
mask(cell2mat(somata)) = false;

trunkAgglos = cellfun( ...
    @(segIds) segIds(mask(segIds)), ...
    dendrites, 'UniformOutput', false);

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
    dendrites, 'UniformOutput', false);

spineCount = cellfun(@numel, spineIds);
spineDensity = spineCount ./ trunkLensUm;

%% Plot spine densities
minSynCount = 10;
maxSpinesPerUm = 0.4;
dendSubMask = ...
    (conn.denMeta.targetClass ~= 'Somata') ...
  & (conn.denMeta.synCount >= minSynCount);

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
    spineDensity(dendSubMask), binEdges, ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);

ax.TickDir = 'out';
xlabel(ax, 'Spine density (µm^{-1})');
ylabel(ax, 'Dendrites');

legend(ax, ...
    'All dendrites', ...
    sprintf('All dendrites with >= %d synapses', minSynCount));

annotation(fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'HorizontalAlignment', 'center', ...
    'EdgeColor', 'none');

%% Export examples to webKNOSSOS
% Very rough threshold based on table 2 from
% Kawaguchi, Karuba, Kubota (2006) Cereb Cortex

randIds = find( ...
    dendSubMask & ...
    spineDensity <= maxSpinesPerUm);
randIds = randIds(randperm(numel(randIds)));
randIds = randIds(1:25);

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

randNodes = cellfun( ...
    @(segIds) segPoints(segIds, :), ...
    dendrites(randIds), 'UniformOutput', false);
randNames = arrayfun( ...
    @(idx, id, rho) sprintf( ...
        '%0*d. Dendrite %g. %.2f spines / µm', ...
        ceil(log10(1 + numel(randIds))), idx, id, rho), ...
    reshape(1:numel(randIds), [], 1), randIds, spineDensity(randIds), ...
    'UniformOutput', false);

skel = Skeleton.fromMST(randNodes, param.raw.voxelSize, skel);
skel.names = randNames;

skel.write(fullfile(outDir, 'smooth-dendrite-candidates.nml'));
