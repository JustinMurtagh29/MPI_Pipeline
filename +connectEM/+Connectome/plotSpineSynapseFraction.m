% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

minSynCount = 10;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

% [conn, syn] = connectEM.Connectome.load(param, connFile);

% HACK(amotta): Switch to "official" connectome loader once the
% intersynpase file becomes available.
conn = load(connFile);
syn = load(conn.info.param.synFile);
conn.axonMeta = connectEM.Axon.completeSynapseMeta(param, conn, syn);

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
