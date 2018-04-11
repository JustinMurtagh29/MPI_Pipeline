% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

% See `+connectEM/addSegmentClassInformation.m`
% Commit hash `4cc2183581ba578585eaacbc0ebc564d4b05aded`
typeEmAggloFile = fullfile(rootDir, 'segmentAggloPredictions.mat.20170523');
typeEmSegFile = fullfile(rootDir, 'segmentPredictions.mat');

% File paths taken from `+connectEM/aggloPreprocessing.m`
% Commit hash `2aac567ca0a704f586ed0e899ea5723a24cd0533`

% Axons from `/gaba/scratch/mberning/aggloGridSearch6/6_01_00046/metricsFinal.mat`
% This in turn is based on `/gaba/scratch/mberning/aggloGridSearch/search05_00564.mat`
% SHA256 hash `a59602e74ea231fb55ca3e1c4b30e5f1da0c7469b34bceed0d16d37133fc917f`
% Extracted `probThresholdAxon` variable
connectEmThreshAxon = 0.97;

% Dendrites from `/gaba/scratch/mberning/aggloGridSearch/search03_00514.mat`
% SHA256 hash `610c6a3f1c09c8f6a7c4708bcd64e6691ea6d6a0ae2fbde84c7d20b3e919e9bf`
% Extracted `probThresholdDendrite` variable
connectEmThreshDend = 0.98;

% See amotta's `+L4/+Spine/+Head/buildAgglos.m`
% Commit hash `36012fb00b88c2d16a3cb4383f32adbdf99370f1`
shThresh = 0.5;

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

%% Plot connectEM / spine head probabilities
binEdges = linspace(0, 1, 21);

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [410, 290];

ax = axes(fig);
hold(ax, 'on');

colors = ax.ColorOrder(1:2, :);

connectEmHist = histogram(ax, ...
    connectEmScores, binEdges, ...
    'DisplayStyle', 'stairs', ...
    'EdgeColor', colors(1, :), ...
    'LineWidth', 2);
typeEmHist = histogram(ax, ...
    typeEmScores(:, 4), binEdges, ...
    'DisplayStyle', 'stairs', ...
    'EdgeColor', colors(2, :), ...
    'LineWidth', 2);

% Plot thresholds
ax.TickDir = 'out';
ax.XLim = binEdges([1, end]);

ax.YScale = 'log';
ax.YLim(1) = exp(-0.1);

ax.YMinorTick = 'off';
ax.YTick = 10 .^ (0:2:8);

plotThresh(ax, connectEmHist, connectEmThreshAxon);
plotThresh(ax, connectEmHist, connectEmThreshDend);
plotThresh(ax, typeEmHist, shThresh);

xlabel(ax, 'Probability');
ylabel(ax, 'Segments / edges');

title(ax, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Utilities
function p = plotThresh(ax, hist, thresh)
    binId = find(hist.BinEdges > thresh, 1);
    binCount = hist.Values(binId - 1);
    
    xVals = repelem(thresh, 1, 2);
    yVals = [ax.YLim(1), binCount];
    
    p = plot( ...
        ax, xVals, yVals, ...
        'Color', hist.EdgeColor, ...
        'LineStyle', '--', ...
        'LineWidth', 2);
end
