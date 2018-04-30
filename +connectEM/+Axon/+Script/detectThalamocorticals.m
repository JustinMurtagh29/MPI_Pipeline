% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

% These axons were manually identified based on the following tracings:
% * TC: https://webknossos.brain.mpg.de/annotations/Explorational/5ae5fff32700006943b3e2a8
% * CC: https://webknossos.brain.mpg.de/annotations/Explorational/5ae6c9c9270000c253b3f60b
tcAxonIds = [ ...
    15197, 4391, 185, 6758, 5761, ...
    4793, 9960, 5712, 7195, 2162];
ccAxonIds = [ ...
    20115, 9052, 19696, 21957, 5040, ...
    13742, 7553, 28847, 20707, 18837];

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

[connDir, connName] = fileparts(connFile);
interSynFile = fullfile(connDir, sprintf('%s_intersynapse.mat', connName));

conn = load(connFile);

syn = load(conn.info.param.synFile);
conn.axonMeta = connectEM.Axon.completeSynapseMeta(param, conn, syn);

interSyn = load(interSynFile);
conn.axonMeta.pathLen = nan(size(conn.axons));
conn.axonMeta.pathLen(interSyn.axonIds) = interSyn.axonPathLens / 1E3;

% Precompute useful quantities
conn.axonMeta.fullPriSpineSynDens = ...
    conn.axonMeta.fullPriSpineSynCount ...
 ./ conn.axonMeta.pathLen;
conn.axonMeta.fullPriSpineSynFrac = ...
    conn.axonMeta.fullPriSpineSynCount ...
 ./ conn.axonMeta.fullSynCount;

%% Plot synapse densities for GT axons
binEdges = linspace(0, 0.4, 21);

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [400, 400];

ax = axes(fig);
hold(ax, 'on');
axis(ax, 'square');

histogram( ...
    ax, conn.axonMeta.fullPriSpineSynDens(tcAxonIds), ...
    'BinEdges', binEdges, 'DisplayStyle', 'stairs', ...
    'LineWidth', 2, 'FaceAlpha', 1);
histogram( ...
    ax, conn.axonMeta.fullPriSpineSynDens(ccAxonIds), ...
    'BinEdges', binEdges, 'DisplayStyle', 'stairs', ...
    'LineWidth', 2, 'FaceAlpha', 1);

ax.TickDir = 'out';
xlim(ax, binEdges([1, end]));
yticks(ax, 0:ceil(ax.YLim(2)));

leg = legend(ax, ...
    'Thalamocortical', ...
    'Corticocortical', ...
    'Location', 'NorthEast');
leg.Box = 'off';

xlabel(ax, 'Spine synapse density (µm^{-1})');
ylabel(ax, 'Axons');

title( ...
    ax, {info.filename, info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Extrapolate to all excitatory axons
thresh = 0.19;
excAxonIds = find( ...
    conn.axonMeta.synCount >= 10 ...
  & conn.axonMeta.fullPriSpineSynFrac >= 0.5);

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [400, 400];

ax = axes(fig);
axis(ax, 'square');
hold(ax, 'on');

histogram( ...
    ax, conn.axonMeta.fullPriSpineSynDens(excAxonIds), ...
    'BinEdges', binEdges, 'DisplayStyle', 'stairs', ...
    'LineWidth', 2, 'FaceAlpha', 1);
plot( ...
    ax, [thresh, thresh], ax.YLim, ...
    'Color', 'black', 'LineStyle', '--');

ax.TickDir = 'out';
xlim(ax, binEdges([1, end]));

xlabel(ax, 'Spine synapse density (µm^{-1})');
ylabel(ax, 'Excitatory axons');

title( ...
    ax, {info.filename, info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
