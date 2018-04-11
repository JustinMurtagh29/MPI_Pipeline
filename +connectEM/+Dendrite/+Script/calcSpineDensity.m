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
spinyFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_spine_attachment.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');

nmlDir = '';
isoDir = '/tmpscratch/amotta/l4/2018-04-11-smooth-dendrite-isosurfaces';

% WKW is faster
segParam = struct;
segParam.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
segParam.backend = 'wkwrap';

% Very rough threshold based on table 2 from
% Kawaguchi, Karuba, Kubota (2006) Cereb Cortex
maxSpinesPerUm = 0.4;
minSynCount = 10;

% Coordinates of interneuron (IN) somata
% KAMIN cell #28 is missing because its soma is not within the dataset and
% is thus not among our set of "whole cells".
inSomaPos = [ ...
    453, 4416, 2507;  % KAMIN cell #21
    3858, 5265, 395]; % KAMIN cell #29
    
info = Util.runInfo();

%% Loading data
Util.log('Loading data');
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
param.seg = segParam;

maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);

% Dendrite prior to spine attachment
trunks = load(trunkFile);
trunks = trunks.dendrites(trunks.indBigDends);
trunks = Agglo.fromSuperAgglo(trunks);

% Spine heads
spiny = load(spinyFile);
spineHeads = spiny.shAgglos;
spinyWholeCells = spiny.dendAgglos(spiny.indWholeCells);

% Connectome (and dendrites after spine attachment)
conn = load(connFile);
dendrites = conn.dendrites;

somata = (conn.denMeta.targetClass == 'Somata');
somata = conn.dendrites(somata);

%% Remove spine heads and necks from dendrites
Util.log('Removing spine heads and necks');
mask = false(maxSegId, 1);
mask(cell2mat(trunks)) = true;
mask(cell2mat(spineHeads)) = false;
mask(cell2mat(somata)) = false;

trunkAgglos = cellfun( ...
    @(segIds) segIds(mask(segIds)), ...
    dendrites, 'UniformOutput', false);

trunkSuperAgglos = cellfun( ...
    @(segIds) segPoints(segIds, :), ...
    trunkAgglos, 'UniformOutput', false);
trunkSuperAgglos = struct('nodes', trunkSuperAgglos);

%% Find IN whole cells
somaPos = cell2mat(cellfun( ...
    @(ids) mean(segPoints(ids, :), 1), ...
    somata, 'UniformOutput', false));

inSomaIds = pdist2( ...
    param.raw.voxelSize .* somaPos, ...
    param.raw.voxelSize .* inSomaPos);
[~, inSomaIds] = min(inSomaIds, [], 1);

inSomaSegIds = somata(inSomaIds);
inSomaLUT = Agglo.buildLUT(maxSegId, inSomaSegIds);

inSpinyWholeCellIds = cellfun( ...
    @(ids) setdiff(inSomaLUT(ids), 0), ...
    spinyWholeCells, 'UniformOutput', false);
inSpinyWholeCellIds = find(~cellfun(@isempty, inSpinyWholeCellIds));
assert(numel(inSpinyWholeCellIds) == numel(inSomaIds));

%% Calculate path lengths
Util.log('Calculating path length');
trunkLensUm = Superagglos.mstLength( ...
    trunkSuperAgglos, param.raw.voxelSize);
trunkLensUm = trunkLensUm / 1E3;

%% Calculate spine density
Util.log('Calculating spine density');
shLUT = Agglo.buildLUT(maxSegId, spineHeads);
spineIds = cellfun( ...
    @(segIds) setdiff(shLUT(segIds), 0), ...
    dendrites, 'UniformOutput', false);

spineCount = cellfun(@numel, spineIds);
spineDensity = spineCount ./ trunkLensUm;

%% Plot spine densities
candMask = ...
    (conn.denMeta.targetClass ~= 'Somata') ...
  & (conn.denMeta.targetClass ~= 'WholeCell') ...
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
    spineDensity(candMask), binEdges, ...
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

%% Select smooth dendrites
smoothIds = find(candMask & (spineDensity <= maxSpinesPerUm));

%% Export examples to webKNOSSOS
if ~isempty(nmlDir)
    randIds = smoothIds;
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

    skel.write(fullfile(nmlDir, 'smooth-dendrite-candidates.nml'));
end

%% Generate isosurfaces
if ~isempty(isoDir)
    Util.log('Generating isosurfaces');
    Visualization.exportAggloToAmira( ...
        dendrites(smoothIds), isoDir, ...
        'smoothSizeHalf', 4, ...
        'smoothWidth', 8, ...
        'reduce', 0.05);
end

Util.log('Done!');
