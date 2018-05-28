% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

targetClasses = { ...
    'Somata', 'SO'; ...
    'ProximalDendrite', 'PD'; ...
    'SmoothDendrite', 'SD'; ...
    'ApicalDendrite', 'AD'; ...
    'AxonInitialSegment', 'AIS'; ...
    'OtherDendrite', 'Other'};

targetLabels = targetClasses(:, 2);
targetClasses = targetClasses(:, 1);

minSynPre = 10;
info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, ~, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);
conn = ...
    connectEM.Connectome.prepareForSpecificityAnalysis(conn);

%% Calculate synapse fractions
classConnectome = ...
    connectEM.Connectome.buildClassConnectome( ...
        conn, 'targetClasses', targetClasses);

classConnectome(conn.axonMeta.synCount < minSynPre, :) = [];
classConnectome = classConnectome ./ sum(classConnectome, 2);

%% Plot results
boxGroups = 1:numel(targetClasses);
boxGroups = repmat(boxGroups, size(classConnectome, 1), 1);

% Add jitter
boxJitter = 0.5 * (rand(size(boxGroups)) - 0.5);
boxJitter = boxGroups + boxJitter;

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [720, 300];

ax = axes(fig);
hold(ax, 'on');

% Fake plots for legend
colors = ax.ColorOrder(1:numel(targetClasses), :);
cellfun(@(c) plot(ax, nan, nan, '.', 'Color', c), num2cell(colors, 2));

scatter( ...
    ax, boxJitter(:), classConnectome(:), ...
    [], colors(boxGroups(:), :), '.');

for curIdx = 1:numel(targetClasses)
    curY = mean(classConnectome(:, curIdx));
    curY = repelem(curY, 1, 2);
    
    curX = curIdx + [-0.45, +0.45];
    curColor = colors(curIdx, :);
    
    plot(ax, ...
        curX, curY, ...
        'Color', curColor, ...
        'LineWidth', 2);
end

ax.Box = 'off';
ax.TickDir = 'out';
ax.YScale = 'log';

yticklabels(ax, arrayfun( ...
    @(f) sprintf('%g', f), ...
    yticks(ax), 'UniformOutput', false));
ax.YAxis.MinorTick = 'off';
ylabel(ax, 'Synapse fraction');

xlim(ax, [0, numel(targetClasses)] + 0.5);
xticks(ax, 1:numel(targetClasses));
xticklabels(ax, targetLabels);

title(ax, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
