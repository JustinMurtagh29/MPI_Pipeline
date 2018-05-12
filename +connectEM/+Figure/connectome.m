% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

minSynCount = 10;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = connectEM.Connectome.load(param, connFile);

%%
axonMeta = conn.axonMeta;
axonMeta(axonMeta.synCount < minSynCount, :) = [];

dendMeta = conn.denMeta;
dendMeta(dendMeta.synCount < minSynCount, :) = [];

% Separate between exc. and inh. cells
inMask = dendMeta.isInterneuron;
somaMask = dendMeta.targetClass == 'Somata';
dendMeta.targetClass(somaMask &  inMask) = 'SomaInh';
dendMeta.targetClass(somaMask & ~inMask) = 'SomaExc';

% Also rename "whole cells" to "proximal dendrites"
proxDendMask = dendMeta.targetClass == 'WholeCell';
dendMeta.targetClass(proxDendMask &  inMask) = 'ProximalDendriteInh';
dendMeta.targetClass(proxDendMask & ~inMask) = 'ProximalDendriteExc';

conn = conn.connectome;
conn(~ismember(conn.edges(:, 1), axonMeta.id), :) = [];
conn(~ismember(conn.edges(:, 2), dendMeta.id), :) = [];

%% Group by classes
axonClasses = { ...
    'Corticocortical'; ...
    'Thalamocortical'; ...
    'Inhibitory';};
dendClasses = { ...
    'SomaExc'; ...
    'SomaInh';
    'ProximalDendriteExc';
    'ProximalDendriteInh';
    'SmoothDendrite'; ...
    'ApicalDendrite'; ...
    'AxonInitialSegment'; ...
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
fig.Color = 'black';
fig.Position(3:4) = 1200;

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

ax.XAxis.Color = 'white';
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

ax.YAxis.Color = 'white';
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = yMinorTicks;

yticks(ax, yTicks);
yticklabels(ax, axonClasses);
ylim(ax, yMinorTicks([1, end]));

ax.Position = [0.15, 0.15, 0.7, 0.7];

annotation( ...
    fig, 'textbox', [0, 0.9, 1, 0.1], 'Color', 'white', ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');
