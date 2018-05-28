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

plotTargetClasses = { ...
    'ApicalDendrite', ...
    'OtherDendrite', ...
    'SmoothDendrite'};

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

%% Calculating data
dendT = conn.denMeta;
dendT(dendT.synCount < minSynPost, :) = [];
dendT = dendT(ismember(dendT.targetClass, plotTargetClasses), :);

dendAgglos = conn.dendrites(dendT.id);

%% Smooth dendrites
Util.log('Calculating spine density');

dendT.trunkLen = ...
    connectEM.Dendrite.calculatePathLengths(param, dendAgglos, trunks);
dendT.spineCount = ...
    connectEM.Dendrite.calculateSpineCount(param, dendAgglos, shAgglos);
dendT.spineDensity = dendT.spineCount ./ (dendT.trunkLen / 1E3);

%% Plot
binEdges = linspace(0, 2, 41);

fig = figure();
fig.Color = 'white';

ax = axes(fig);
axis(ax, 'square');
hold(ax, 'on');

plotHist = @(data) ...
    histogram( ...
        ax, data, ...
        'BinEdges', binEdges, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2, ...
        'FaceAlpha', 1);
    
plotHist(dendT.spineDensity);
plotHist(dendT.spineDensity(dendT.targetClass == 'ApicalDendrite', :));
plotHist(dendT.spineDensity(dendT.targetClass == 'SmoothDendrite', :));

plot( ...
    ax, repelem(maxSpinesPerUm, 1, 2), ax.YLim, ...
    'Color', 'black', 'LineStyle', '--');

ax.TickDir = 'out';
xlabel(ax, 'Spine density (µm^{-1})');
ylabel(ax, 'Dendrites');

leg = { ...
    'All dendrites', ...
    'Apical dendrites', ...
    'Smooth dendrites'};
leg = legend(ax, leg, 'Location', 'NorthEast');
leg.Box = 'off';

title( ...
    ax, {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
