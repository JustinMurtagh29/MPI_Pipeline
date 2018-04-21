% This script attempts to calculate the dendrite path length along the
% trunk. This is done by calculating the MST after removing spine head and
% neck segments.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
trunkFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v2.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v2_auto.mat');
dendFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_03.mat');
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');

% Very rough threshold based on table 2 from
% Kawaguchi, Karuba, Kubota (2006) Cereb Cortex
maxSpinesPerUm = 0.4;
minSynCount = 10;
    
info = Util.runInfo();

%% Loading data
Util.log('Loading data');
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

% Dendrite prior to spine attachment
trunks = load(trunkFile);
trunks = trunks.dendrites(trunks.indBigDends);
trunks = Agglo.fromSuperAgglo(trunks);

% Spine heads
shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

% Dendrites + somata
dend = load(dendFile);

somaAgglos = dend.dendrites(dend.indSomata);
somaAgglos = Agglo.fromSuperAgglo(somaAgglos);

dendIds = dend.indAIS | dend.indSomata | dend.indWholeCells;
dendIds = find(dend.indBigDends & ~dendIds);

dend = dend.dendrites(dendIds);
dend = Agglo.fromSuperAgglo(dend);

% Synapses
syn = load(synFile);
syn = syn.synapses;

%% Calculate spine density
Util.log('Calculating spine density');

trunkLensUm = ...
    connectEM.Dendrite.calculatePathLengths(param, dend, trunks);
trunkLensUm = trunkLensUm / 1E3;

spineDensity = ...
    connectEM.Dendrite.calculateSpineDensity( ...
        param, dend, trunkLensUm, shAgglos);
    
%% Determine number of synapses per dendrite
maxSegId = Seg.Global.getMaxSegId(param);
dendLUT = Agglo.buildLUT(maxSegId, dend);

syn.dendIds = cellfun( ...
    @(ids) reshape(setdiff(dendLUT(ids), 0), [], 1), ...
    syn.postsynId, 'UniformOutput', false);
dendSynCount = accumarray( ...
    cell2mat(syn.dendIds), 1, size(dend));

%% Plot spine densities
candMask = dendSynCount >= minSynCount;

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
