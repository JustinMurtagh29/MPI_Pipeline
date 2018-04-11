% This script attempts to calculate the dendrite path length along the
% trunk. This is done by calculating the MST after removing spine head and
% neck segments.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
trunkFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_v2.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_spine_attachment.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');

% Very rough threshold based on table 2 from
% Kawaguchi, Karuba, Kubota (2006) Cereb Cortex
maxSpinesPerUm = 0.4;
minSynCount = 10;
    
info = Util.runInfo();

%% Loading data
Util.log('Loading data');
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);

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

somaAgglos = conn.denMeta.targetClass == 'Somata';
somaAgglos = conn.dendrites(somaAgglos);

%% Calculate spine density
Util.log('Calculating spine density');

trunkLensUm = ...
    connectEM.Dendrite.calculatePathLengths( ...
        param, dendrites, trunks, shAgglos, somaAgglos);
trunkLensUm = trunkLensUm / 1E3;

spineCount = ...
    Agglo.buildLUT(maxSegId, shAgglos);
spineCount = cellfun( ...
    @(segIds) setdiff(spineCount(segIds), 0), ...
    dendrites, 'UniformOutput', false);
spineCount = cellfun(@numel, spineCount);

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
