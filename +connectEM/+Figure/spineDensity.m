% Based on
%   +connectEM/+Dendrite/+Script/classify.m
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
trunkFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

% Very rough threshold based on table 2 from
% Kawaguchi, Karuba, Kubota (2006) Cereb Cortex
maxSpinesPerUm = 0.4;
minSynPost = 10;

info = Util.runInfo();

%% Loading data
Util.log('Loading data');
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = connectEM.Connectome.load(param, connFile);

% Dendrite trunks (i.e., prior to spine attachment)
trunks = load(trunkFile);
trunks = trunks.dendrites(trunks.indBigDends);
trunks = Agglo.fromSuperAgglo(trunks);

% Spine heads
shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

dendMask = {'AxonInitialSegment', 'Somata', 'WholeCell'};
dendMask = ismember(conn.denMeta.targetClass, dendMask);
dendMask = (conn.denMeta.synCount >= minSynPost) & ~dendMask;
dendAgglos = conn.dendrites(dendMask);

%% Smooth dendrites
Util.log('Calculating spine density');

trunkLens = ...
    connectEM.Dendrite.calculatePathLengths(param, dendAgglos, trunks);
spineCounts = ...
    connectEM.Dendrite.calculateSpineCount(param, dendAgglos, shAgglos);
spineDens = spineCounts ./ (trunkLens / 1E3);

%% Plot
binEdges = linspace(0, 2, 41);

fig = figure();
fig.Color = 'white';

ax = axes(fig);
axis(ax, 'square');
hold(ax, 'on');

histogram( ...
    ax, spineDens, ...
    'BinEdges', binEdges, ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);
plot( ...
    ax, repelem(maxSpinesPerUm, 1, 2), ax.YLim, ...
    'Color', 'black', 'LineStyle', '--');

ax.TickDir = 'out';
xlabel(ax, 'Spine density (µm^{-1})');
ylabel(ax, 'Dendrites');

title( ...
    ax, {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
