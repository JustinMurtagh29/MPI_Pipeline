% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

% See `+connectEM/addSegmentClassInformation.m`
% Commit hash `4cc2183581ba578585eaacbc0ebc564d4b05aded`
typeEmAggloFile = fullfile(rootDir, 'segmentAggloPredictions.mat.20170523');
typeEmSegFile = fullfile(rootDir, 'segmentPredictions.mat');

% As in Benedikt's +L4/updateParamsToNewestFiles.m
% Commit hash `590d8538d65463151d43491e2446e25ca11dd5f6`
synEmGraphFile = fullfile(rootDir, 'graphNew.mat');
synEmScoreFile = fullfile(rootDir, 'globalSynScores.mat');

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

% Loading SynEM scores
param.svg.graphFile = synEmGraphFile;
param.svg.synScoreFile = synEmScoreFile;
graph = Seg.IO.loadGraph(param, false);

synEmScores = graph.synScores;
synEmScores(isnan(graph.borderIdx), :) = [];
clear graph;

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

%% Plot SynEM scores
binEdges = linspace(-20, 5, 51);

rng(0);
randSynEmScores = randperm(size(synEmScores, 1));
randSynEmScores = randSynEmScores(1:1e4);
randSynEmScores = synEmScores(randSynEmScores, :);

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [680, 680];

ax = axes(fig);
ax.Position = [0.25, 0.25, 0.65, 0.65];

scatter(ax, ...
    randSynEmScores(:, 1), ...
    randSynEmScores(:, 2), '.');
axis(ax, 'square');

xticks(ax, []);
xlim(ax, binEdges([1, end]));
yticks(ax, []);
ylim(ax, binEdges([1, end]));

% Bottom
axOne = axes(fig);
axOne.Position = [0.25, 0.1, 0.65, 0.1];

histogram(axOne, ...
    randSynEmScores(:, 1), binEdges, ...
    'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);

axOne.TickDir = 'out';
axOne.YDir = 'reverse';
xlim(axOne, binEdges([1, end]));
xlabel('SynEM score 1');

% Left
axTwo = axes(fig);
axTwo.Position = [0.1, 0.25, 0.1, 0.65];

histogram(axTwo, ...
    randSynEmScores(:, 2), binEdges, ...
    'Normalization', 'probability', ...
    'Orientation', 'horizontal', ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);

axTwo.TickDir = 'out';
axTwo.XDir = 'reverse';
ylim(axTwo, binEdges([1, end]));
ylabel('SynEM score 2');

% Make sure both histogram have same ranges
axOne.YLim(2) = max(axOne.YLim(2), axTwo.XLim(2));
axTwo.XLim(2) = max(axOne.YLim(2), axTwo.XLim(2));

% Annotation
annotation(fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'String', {info.filename; info.git_repos{1}.hash});
