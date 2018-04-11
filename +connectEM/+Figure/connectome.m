% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');

minSynCount = 10;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[~, connName] = fileparts(connFile);
conn = connectEM.Connectome.load(param, connName);
connBack = conn;

%%
axonMeta = conn.axonMeta;
axonMeta(axonMeta.synCount < minSynCount, :) = [];

dendMeta = conn.denMeta;
dendMeta(dendMeta.synCount < minSynCount, :) = [];

conn = conn.connectome;
conn(~ismember(conn.edges(:, 1), axonMeta.id), :) = [];
conn(~ismember(conn.edges(:, 2), dendMeta.id), :) = [];

%% Group by classes
axonClasses = { ...
    'Corticocortical'; ...
    'Inhibitory'; ...
    'Thalamocortical'; ...
    'Other'};
dendClasses = { ...
    'Somata'; ...
    'WholeCell'; ...
    'AxonInitialSegment'; ...
    'ApicalDendrite'; ...
    'SmoothDendrite'; ...
    'OtherDendrite'};

rng(0);
randIds = randperm(numel(axonMeta.id));
axonMeta = sortrows(axonMeta, 'id');
axonMeta = axonMeta(randIds, :);

[~, axonMeta.classIdx] = ismember( ...
    axonMeta.axonClass, axonClasses);

assert(all(axonMeta.classIdx));
axonMeta = sortrows(axonMeta, 'classIdx');

rng(0);
randIds = randperm(numel(dendMeta.id));
dendMeta = sortrows(dendMeta, 'id');
dendMeta = dendMeta(randIds, :);

[~, dendMeta.classIdx] = ismember( ...
    dendMeta.targetClass, dendClasses);

assert(all(dendMeta.classIdx));
dendMeta = sortrows(dendMeta, 'classIdx');

%% Plot
colorMap = parula(5);

conn.synCount = cellfun(@numel, conn.synIdx);
conn.coord = nan(size(conn.edges));
[~, conn.coord(:, 1)] = ismember(conn.edges(:, 1), axonMeta.id);
[~, conn.coord(:, 2)] = ismember(conn.edges(:, 2), dendMeta.id);
conn = sortrows(conn, 'synCount', 'ascend');

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = 1100;

ax = axes(fig);
h = scatter(ax, conn.coord(:, 2), conn.coord(:, 1), '.');
h.CData = colorMap(min(conn.synCount, size(colorMap, 1)), :);

ax.TickDir = 'out';
ax.YDir = 'reverse';
ax.XAxisLocation = 'top';

ax.Color = 'black';

% Target classes
xMinorTicks = accumarray( ...
    dendMeta.classIdx, ...
    1:numel(dendMeta.id), ...
    size(dendClasses), @max);
xMinorTicks = [1; xMinorTicks(:)];

xTicks = ...
    xMinorTicks(1:(end - 1)) ...
  + xMinorTicks(2:end);
xTicks = xTicks / 2;

ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = xMinorTicks;

xticks(ax, xTicks);
xticklabels(ax, dendClasses);
xtickangle(ax, 90);
xlim(ax, xMinorTicks([1, end]));

% Axon classes
yMinorTicks = accumarray( ...
    axonMeta.classIdx, ...
    1:numel(axonMeta.id), ...
    size(axonClasses), @max);
yMinorTicks = [1; yMinorTicks(:)];

yTicks = ...
    yMinorTicks(1:(end - 1)) ...
  + yMinorTicks(2:end);
yTicks = yTicks / 2;

ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = yMinorTicks;

yticks(ax, yTicks);
yticklabels(ax, axonClasses);
ylim(ax, yMinorTicks([1, end]));

ax.Position = [0.15, 0.15, 0.7, 0.7];

% Histogram over incoming synapses
axDend = axes(fig);
axDend.Position = ax.Position;
axDend.Position(2:2:end) = [0.05, ax.Position(2) - 0.05];

dendSynCount = accumarray( ...
    conn.coord(:, 2), conn.synCount, size(dendMeta.id));
histogram(axDend, ...
    'BinCount', dendSynCount, ...
    'BinEdges', 0:1:numel(dendSynCount), ...
    'EdgeColor', ax.ColorOrder(1, :));

axDend.YScale = 'log';
axDend.YDir = 'reverse';
axDend.TickDir = 'out';
axDend.YMinorTick = 'off';

xlim(axDend, [0, numel(dendSynCount)]);
xticks(axDend, [xticks(axDend), axDend.XLim(2)]);
ylabel(axDend, 'Synapses');

% Histogram over outgoing synapses
axAxon = axes(fig);
axAxon.Position = ax.Position;
axAxon.Position(1) = sum(ax.Position(1:2:end));
axAxon.Position(3) = 0.95 - axAxon.Position(1);

axonSynCount = accumarray( ...
    conn.coord(:, 1), conn.synCount, size(axonMeta.id));
histogram(axAxon, ...
    'BinCounts', axonSynCount, ...
    'BinEdges', 0:1:numel(axonSynCount), ...
    'Orientation', 'horizontal', ...
    'EdgeColor', ax.ColorOrder(1, :));

axAxon.TickDir = 'out';
axAxon.XScale = 'log';
axAxon.XAxisLocation = 'top';
axAxon.XMinorTick = 'off';
axAxon.YDir = 'reverse';
axAxon.YAxisLocation = 'right';

ylim(axAxon, [0, numel(axonSynCount)]);
yticks(axAxon, [yticks(axAxon), numel(axonSynCount)]);
xlabel('Synapses');

annotation( ...
    fig, 'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');
