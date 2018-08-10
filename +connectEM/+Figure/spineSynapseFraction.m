% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

minSynCount = 10;
threshLines = [0.2, 0.5];

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = connectEM.Connectome.load(param, connFile);

%% Prepare axon meta data
axonMeta = conn.axonMeta;

% Remove axons with too few synapses
axonMeta(axonMeta.synCount < minSynCount, :) = [];

axonMeta.fullPriSpineSynFrac = ...
    axonMeta.fullPriSpineSynCount ...
 ./ axonMeta.fullSynCount;

%% Plot primary spine synapse fraction
binEdges = linspace(0, 1, 21);

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [290, 210];

ax = axes(fig);
axis(ax, 'square')
hold(ax, 'on');

histogram(ax, ...
    axonMeta.fullPriSpineSynFrac, binEdges, ...
    'DisplayStyle', 'stairs', 'LineWidth', 2, 'FaceAlpha', 1);

lines = arrayfun(@(x) plot(ax, [x, x], ax.YLim), threshLines);
set(lines, 'Color', 'black', 'LineWidth', 2', 'LineStyle', '--');

ax.TickDir = 'out';
ax.XLim = binEdges([1, end]);

xlabel(ax, 'Spine synapse fraction');
ylabel(ax, 'Axons');

title(ax, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
