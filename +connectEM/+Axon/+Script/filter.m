% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_18_b.mat');

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segmentMeta = struct;
segmentMeta.maxSegId = Seg.Global.getMaxSegId(param);
segmentMeta = connectEM.addSegmentClassInformation(param, segmentMeta);

axon = load(axonFile);

%% Calculate per-axon glia score
axonIds = find(axon.indBigAxons);

axons = axon.axons(axonIds);
axons = Agglo.fromSuperAgglo(axons);

gliaScore = segmentMeta.gliaProb;
gliaScore = cellfun(@(ids) median(gliaScore(ids), 'omitnan'), axons);

%% Histogram of glia score
binEdges = linspace(0, 1, 51);

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [460, 430];

ax = axes(fig);

histogram(ax, ...
    gliaScore, binEdges, ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);

ax.Box = 'off';
ax.TickDir = 'out';
ax.YScale = 'log';

ax.XLim = binEdges([1, end]);
ax.YLim(1) = 0;

axis(ax, 'square');
xlabel(ax, 'Median glia probability');
ylabel(ax, 'Axons');

title(ax, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
