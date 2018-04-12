% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered_classified.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');

minSynCount = 10;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn] = ...
    connectEM.Connectome.load(param, connFile, synFile);

%% Prepare axon meta data% Remove axons with too few synapses
axonMeta = conn.axonMeta;
axonMeta(axonMeta.synCount < minSynCount, :) = [];

axonMeta.fullPriSpineSynFrac = ...
    axonMeta.fullPriSpineSynCount ...
 ./ axonMeta.fullSynCount;

%% Plot primary spine synapse fraction
binEdges = linspace(0, 1, 21);

fig = figure();
fig.Color = 'white';

ax = axes(fig);
axis(ax, 'square')
hold(ax, 'on');

histogram(ax, ...
    axonMeta.fullPriSpineSynFrac, binEdges, ...
    'DisplayStyle', 'stairs', 'LineWidth', 2, 'FaceAlpha', 1);

ax.TickDir = 'out';
xlim(ax, binEdges([1, end]));
xlabel(ax, 'Spine synapse fraction');
ylabel(ax, 'Axons');

title(ax, ...
   {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
