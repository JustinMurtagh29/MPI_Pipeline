% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

targetClasses = { ...
    'Somata', 'SO'; ...
    'ProximalDendrite', 'PD'; ...
    'ApicalDendrite', 'AD'; ...
    'SmoothDendrite', 'SD'; ...
    'AxonInitialSegment', 'AIS'};
targetLabels = targetClasses(:, 2);
targetClasses = targetClasses(:, 1);

minSynPre = 10;

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, ~, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);
[conn, axonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, axonClasses, 'minSynPre', minSynPre);

%% Prepare data
axonClasses = axonClasses(1:3);

axonClasses(1).tag = 'Exc';
axonClasses(2).tag = 'Inh';
axonClasses(3).tag = 'TC';

axonClasses = ...
    connectEM.Connectome.buildAxonSpecificityClasses(conn, axonClasses);
[fracSpec, fracOfSpecOntoTarget, fracOfSpecOntoTargetBars] = ...
    connectEM.Specificity.calcAxonFractions(axonClasses, targetClasses);

disp(array2table( ...
    horzcat(fracSpec, fracOfSpecOntoTarget), ...
    'VariableNames', cat(1, 'Overall', targetClasses), ...
    'RowNames', {axonClasses.tag}));

%% Plot fraction of axons specific
clear cur*;
plotClasses = 1:numel(axonClasses);

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [150, 161];

ax = axes(fig);
bar(ax, ...
    1:numel(plotClasses), ...
    fracSpec(plotClasses), ...
    'EdgeColor', 'black', 'FaceColor', 'black');

ax.XLim = [0.5, numel(plotClasses) + 0.5];
ax.YLim = [0, 1];

ax.Box = 'off';
ax.TickDir = 'out';
ax.XTick = 1:numel(plotClasses);
ax.XTickLabel = {axonClasses(plotClasses).tag};

ylabel(ax, {'Fraction of'; 'axons specific'});

title(ax, { ...
    info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Plot distribution of specificities over target classes
clear cur*;

curBarWidth = 0.8;
curHistWidth = 0.09;
plotClasses = 1:2;

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [280, 280];

ax = axes(fig);
hold(ax, 'on');
ax.YAxisLocation = 'right';

plotData = fracOfSpecOntoTargetBars(plotClasses, :);
plotData = flip(plotData, 2);

allBars = bar(ax, plotData, 'stacked', 'BarWidth', curBarWidth);

colors = flip(ax.ColorOrder(1:numel(targetClasses), :), 1);
mixedColors = (colors(1:(end - 1), :) .^ 2 + colors(2:end, :) .^ 2) / 2;
mixedColors = cat(1, sqrt(mixedColors), zeros(1, 3));

colors = cat(1, transpose(colors), transpose(mixedColors));
colors = transpose(reshape(colors, 3, []));
colors = num2cell(colors(1:(end - 1), :), 2);

[allBars.FaceColor] = deal(colors{:});
[allBars.EdgeColor] = deal('none');

excSpecs = axonClasses(1).specs;
excOff = fracOfSpecOntoTargetBars(1, :);
excOff = 1 - cumsum(cat(1, 0, excOff(:)));
for curClassId = 1:numel(targetClasses)
    curClassName = targetClasses{curClassId};
    curTcProbs = [];
    
    if isfield(excSpecs, curClassName)
        curTcProbs = excSpecs.(curClassName).axonIds;
        curTcProbs = conn.axonMeta.tcProb(curTcProbs);
        curTcProbs = sort(curTcProbs, 'descend');
    end
    
    curOffX = 1 - curBarWidth / 2;
    curX = (curTcProbs(:) >= 0.6) * curHistWidth;
    curOffY = excOff(1 + 2 * (curClassId - 1) + [0, 1]);
    curY = linspace(curOffY(1), curOffY(2), numel(curX));
    
    curX = curOffX - cat(1, 0, curX(:), 0);
    curY = cat(1, curOffY(1), curY(:), curOffY(2));
    
    plot(ax, curX, curY, 'Color', 'black');
end

xlim(ax, [0.5, numel(plotClasses) + 0.5]);
xticklabels(ax, {axonClasses(plotClasses).tag});

leg = legend( ...
    flip(allBars(1:2:end)), targetLabels, ...
    'Location', 'EastOutside');
leg.Box = 'off';

ax.Box = 'off';
ax.TickDir = 'out';

title( ...
    ax, {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
