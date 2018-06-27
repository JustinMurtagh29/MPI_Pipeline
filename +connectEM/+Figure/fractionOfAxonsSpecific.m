% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

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

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, ~, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

%% Preparing data
axonClasses(1).tag = 'Exc';
axonClasses(2).tag = 'Inh';
axonClasses(3).tag = 'TC';
axonClasses(4).tag = 'CC';

[conn, axonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, axonClasses, 'minSynPre', minSynPre);
axonClasses = ...
    connectEM.Connectome.buildAxonSpecificityClasses(conn, axonClasses);

%% Prepare data
fractionSpecific = nan(size(axonClasses));
fractionOfSpecificOntoTarget = nan( ...
    numel(axonClasses), 2 * numel(targetClasses) - 1);

for curId = 1:numel(axonClasses)
    curAxonClass = axonClasses(curId);
    curAxonIds = curAxonClass.axonIds;
    curSpecs = curAxonClass.specs;
    
    curSpecClasses = fieldnames(curSpecs);
   [~, curIds] = ismember(curSpecClasses, targetClasses);
    assert(all(curIds));
    
    curSpecAxonIds = cellfun( ...
        @(name) curSpecs.(name).axonIds, ...
        curSpecClasses, 'UniformOutput', false);
    
    % Determine fraction of axons with specificity
    curFractionSpecific = unique(cell2mat(curSpecAxonIds));
    curFractionSpecific = numel(curFractionSpecific) / numel(curAxonIds);
    fractionSpecific(curId) = curFractionSpecific;
    
    % Fraction of specific axons per target
   [curA, curB] = ndgrid(1:numel(curIds), 1:numel(curIds));
    curSpecMat = zeros(numel(targetClasses));
    curSpecMat(curIds, curIds) = cellfun( ...
        @(idsOne, idsTwo) numel(intersect(idsOne, idsTwo)), ...
        curSpecAxonIds(curA), curSpecAxonIds(curB));
    
    % Prepare stacked bars
    curDiag = diag(curSpecMat, 0);
    curOff = diag(curSpecMat, 1);
    
    curDiag(1:(end - 1)) = curDiag(1:(end - 1)) - curOff;
    curDiag(2:end) = curDiag(2:end) - curOff;
    
    curDiag = reshape(curDiag, 1, []);
    curOff = reshape(curOff, 1, []);
    
    curProbs = cat(1, curDiag, cat(2, curOff, 0));
    curProbs = curProbs(1:(end - 1));
    curProbs = curProbs / sum(curProbs);
    
    fractionOfSpecificOntoTarget(curId, :) = curProbs;
end

%% Plot fraction of axons specific
plotClasses = [4, 3, 2];

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [150, 161];

ax = axes(fig);
bar(ax, ...
    1:numel(plotClasses), ...
    fractionOfAxonsSpecific(plotClasses), ...
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
plotClasses = [4, 3, 2];

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [190, 160];

ax = axes(fig);
plotData = fractionOfSpecificOntoTarget(plotClasses, :);
plotData = flip(plotData, 2);

allBars = bar(ax, plotData, 'stacked');

colors = flip(ax.ColorOrder(1:numel(targetClasses), :), 1);
mixedColors = (colors(1:(end - 1), :) .^ 2 + colors(2:end, :) .^ 2) / 2;
mixedColors = cat(1, sqrt(mixedColors), zeros(1, 3));

colors = cat(1, transpose(colors), transpose(mixedColors));
colors = transpose(reshape(colors, 3, []));
colors = num2cell(colors(1:(end - 1), :), 2);

[allBars.FaceColor] = deal(colors{:});
[allBars.EdgeColor] = deal('none');

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
