% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

minSynCount = 10;

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = connectEM.Connectome.load(param, connFile);

%% Filter neurites
% Only plot neurites with at least `minSynCount` synapses. Notable
% exception are AIS, which we plot irrespectively of their synapse number.
axonMeta = conn.axonMeta;
axonMeta(axonMeta.synCount < minSynCount, :) = [];

dendMeta = conn.denMeta;
dendMeta( ...
    dendMeta.targetClass ~= 'AxonInitialSegment' ...
  & dendMeta.synCount < minSynCount, :) = [];

% Rename "Somata" to "Soma"
somaMask = dendMeta.targetClass == 'Somata';
dendMeta.targetClass(somaMask) = 'Soma';

% Rename "WholeCell" to "ProximalDendrite"
proxDendMask = dendMeta.targetClass == 'WholeCell';
dendMeta.targetClass(proxDendMask) = 'ProximalDendrite';

conn = conn.connectome;
conn(~ismember(conn.edges(:, 1), axonMeta.id), :) = [];
conn(~ismember(conn.edges(:, 2), dendMeta.id), :) = [];

%% Numbers
numberOfConnections = height(conn) %#ok

%% Group by classes
axonClasses = { ...
    'Corticocortical'; ...
    'Thalamocortical'; ...
    'Inhibitory';};
dendClasses = { ...
    'Soma'; ...
    'ProximalDendrite';
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
fig.Color = 'white';
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
    'EdgeColor', ax.ColorOrder(1, :), ...
    'FaceAlpha', 1);

axDend.YScale = 'log';
axDend.YDir = 'reverse';
axDend.TickDir = 'out';
axDend.YMinorTick = 'off';

axDend.YLim(1) = minSynCount;
axDend.YLim(2) = 10 ^ ceil(log10(max(dendSynCount)));
yTicks = log10(axDend.YLim);
yTicks = 10 .^ (yTicks(1):yTicks(2));

yticks(axDend, yTicks);
yticklabels(axDend, arrayfun( ...
    @(d) sprintf('%d', d), ...
    yTicks, 'UniformOutput', false));

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
    'EdgeColor', ax.ColorOrder(1, :), ...
    'FaceAlpha', 1);

axAxon.TickDir = 'out';
axAxon.XScale = 'log';
axAxon.XAxisLocation = 'top';
axAxon.XMinorTick = 'off';
axAxon.YDir = 'reverse';
axAxon.YAxisLocation = 'right';

axAxon.XLim(1) = minSynCount;
axAxon.XLim(2) = 10 ^ ceil(log10(max(axonSynCount)));
xTicks = log10(axAxon.XLim);
xTicks = 10 .^ (xTicks(1):xTicks(2));

xticks(axAxon, xTicks);
xticklabels(axAxon, arrayfun( ...
    @(d) sprintf('%d', d), ...
    xTicks, 'UniformOutput', false));

ylim(axAxon, [0, numel(axonSynCount)]);
yticks(axAxon, [yticks(axAxon), numel(axonSynCount)]);
xlabel('Synapses');

% Legend
axLeg = axes(fig);
hold(axLeg, 'on');

axLeg.Position = [ ...
    axAxon.Position(1), axDend.Position(2), ...
    axAxon.Position(3), axDend.Position(4)];

% Invisible plot for legend
legPlots = cellfun( ...
    @(color) plot( ...
        axLeg, nan, nan, '.', ...
        'MarkerSize', 12, ...
        'MarkerEdgeColor', color, ...
        'MarkerFaceColor', color), ...
    num2cell(colorMap, 2));

leg = arrayfun( ...
    @num2str, 1:size(colorMap, 1), ...
    'UniformOutput', false);
leg(end) = strcat(leg(end), '+');

leg = legend(axLeg, leg, 'Location', 'SouthEast');
leg.Position = axLeg.Position;
leg.Box = 'off';

axLeg.XAxis.Visible = 'off';
axLeg.YAxis.Visible = 'off';

annotation( ...
    fig, 'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');
