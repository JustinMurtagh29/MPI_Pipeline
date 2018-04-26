% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-classified_spine-syn-clust.mat');

% This corresponds to the axon shown in panel 3c. This was axon 630 in
% commectome `connectome_axons_18_a_ax_spine_syn_clust`.
axonId = 628;
minSynCount = 10;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, ~, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);
[classConnectome, targetClasses] = ...
    connectEM.Connectome.buildClassConnectome(conn);

%% Prepare data
synCount = sum(classConnectome(axonId, :));
somaSynCount = classConnectome(axonId, targetClasses == 'Somata');

nullAxonIds = find(sum(classConnectome, 2) >= minSynCount);
somaProb = sum(classConnectome(nullAxonIds, :), 1);
somaProb = somaProb(targetClasses == 'Somata') / sum(somaProb);

nullProbs = binopdf(0:synCount, synCount, somaProb);

%% Figure
binEdges = 0:(synCount + 1);

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [400, 250];

ax = axes(fig);
hold(ax, 'on');

histogram(ax, ...
    'BinEdges', binEdges, 'BinCounts', nullProbs, ...
    'DisplayStyle', 'stairs', 'LineWidth', 2, 'FaceAlpha', 1);
plot(ax, ...
    0.5 + repelem(somaSynCount, 1, 2), [0, 1], ...
    'Color', 'black', 'LineStyle', '--');

ax.TickDir = 'out';

xlim(ax, binEdges([1, end]));
ylim(ax, [0, 1]);

xlabel(ax, 'Soma synapses');
ylabel(ax, 'p');

xTicks = union(xticks(ax), [somaSynCount, synCount]);
xticklabels(ax, arrayfun( ...
    @num2str, xTicks, 'UniformOutput', false));
xticks(ax, 0.5 + xTicks);
yticks(ax, [0, 0.5, 1]);

title(ax, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
