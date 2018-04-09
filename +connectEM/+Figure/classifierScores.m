% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

% See `+connectEM/addSegmentClassInformation.m`
% Commit hash `4cc2183581ba578585eaacbc0ebc564d4b05aded`
typeEmAggloFile = fullfile(rootDir, 'segmentAggloPredictions.mat.20170523');
typeEmSegFile = fullfile(rootDir, 'segmentPredictions.mat');

% Plot histogram with log Y axis
plotLog = false;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);

% Loading connectEM scores
graph = Graph.load(rootDir);
connectEmScores = graph.prob;
connectEmScores(~graph.borderIdx) = [];
clear graph;

% Loading TypeEM scores
typeEmScores = nan(maxSegId, 4);

typeEmAgglo = load(typeEmAggloFile, '-mat');
colIds = {'axon', 'dendrite', 'astrocyte'};
[~, colIds] = ismember(colIds, typeEmAgglo.classes);
typeEmScores(typeEmAgglo.segId, 1:3) = typeEmAgglo.probs(:, colIds);
clear typeEmAgglo;

typeEmSeg = load(typeEmSegFile, '-mat');
[~, colIds] = ismember({'spinehead'}, typeEmSeg.class);
typeEmScores(typeEmSeg.segId, 4) = typeEmSeg.probs(:, colIds);
clear typeEmSeg;

typeEmScores(any(isnan(typeEmScores), 2), :) = [];

%% Plot all TypeEM probabilities
rng(0);
randTypeEmScores = randperm(size(typeEmScores, 1));
randTypeEmScores = randTypeEmScores(1:1E4);
randTypeEmScores = typeEmScores(randTypeEmScores, :);
randTypeEmScores = sortrows(randTypeEmScores, 4, 'descend');

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [500, 500];

ax = axes(fig);
h = scatter3(ax, ...
    randTypeEmScores(:, 1), ...
    randTypeEmScores(:, 2), ...
    randTypeEmScores(:, 3), '.');

sizeFunc = @(p) (1 + 4 .* p) * 36;
h.SizeData = sizeFunc(randTypeEmScores(:, 4));

colors = ax.ColorOrder(1:2, :);
colorFunc = @(p) (1 - p) .* colors(1, :) + p .* colors(2, :);
h.CData = colorFunc(randTypeEmScores(:, 4));

xlabel(ax, 'Axon');
ylabel(ax, 'Dendrite');
zlabel(ax, 'Astrocyte');

grid(ax, 'off');
view(ax, -45, 27);
ax.TickDir = 'out';

axis(ax, 'equal');
xlim(ax, [0, 1]);
xticks(ax, [0, 0.5, 1]);
ylim(ax, [0, 1]);
yticks(ax, [0, 0.5, 1]);
ylim(ax, [0, 1]);
zticks(ax, [0, 0.5, 1]);

title(ax, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Plot spine head probability histogram
binEdges = linspace(0, 1, 21);

fig = figure();
fig.Color = 'white';
ax = axes(fig);

histogram(ax, ...
    typeEmScores(:, 4), binEdges, ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);

ax.TickDir = 'out';
ax.XLim = binEdges([1, end]);

if plotLog
    ax.YScale = 'log';
    ax.YLim(1) = exp(-0.1);
else
    ax.YLim(1) = 0;
end

xlabel(ax, 'Spine head probability');
ylabel(ax, 'Segments');

title(ax, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Plot connectEM probability histogram
binEdges = linspace(0, 1, 21);

fig = figure();
fig.Color = 'white';
ax = axes(fig);

histogram(ax, ...
    connectEmScores, binEdges, ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);

ax.TickDir = 'out';
ax.XLim = binEdges([1, end]);

if plotLog
    ax.YScale = 'log';
    ax.YLim(1) = exp(-0.1);
else
    ax.YLim(1) = 0;
end

xlabel(ax, 'connectEM probability');
ylabel(ax, 'Edges');

title(ax, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

