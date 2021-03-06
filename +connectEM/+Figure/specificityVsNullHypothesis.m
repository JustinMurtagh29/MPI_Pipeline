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
nullPVals = fliplr(cumsum(fliplr(nullProbs)));

%% Figure
binEdges = (0:(synCount + 1)) - 0.5;

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [425, 275];

ax = axes(fig);
yyaxis(ax, 'left');
hold(ax, 'on');

histogram(ax, ...
    'BinEdges', binEdges, 'BinCounts', nullProbs, ...
    'DisplayStyle', 'stairs', 'LineWidth', 2, 'FaceAlpha', 1);
plot(ax, ...
    repelem(somaSynCount, 1, 2), [0, 1], ...
    'Color', 'black', 'LineStyle', '--');

ax.TickDir = 'out';

ylim(ax, [0, 1]);

ylabel(ax, 'p');
yticks(ax, [0, 0.5, 1]);

yyaxis(ax, 'right');
plot(ax, (0:synCount), nullPVals, 'LineWidth', 2);
ax.YAxis(2).Scale = 'log';

xTicks = union(xticks(ax), [somaSynCount, synCount]);
xticklabels(ax, arrayfun( ...
    @num2str, xTicks, 'UniformOutput', false));
xticks(ax, xTicks);

xlim(ax, binEdges([1, end]));
xlabel(ax, 'Soma synapses');

title(ax, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
