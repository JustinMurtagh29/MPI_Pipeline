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

%% Plotting
plotClasses = [4, 3, 2];
plotData = nan(numel(plotClasses), 2 * numel(targetClasses) + 1);

for curIdx = 1:numel(plotClasses)
    curId = plotClasses(curIdx);
    curAxonClass = axonClasses(curId);
    curAxonIds = curAxonClass.axonIds;
    curSpecs = curAxonClass.specs;
    
    curSpecClasses = fieldnames(curSpecs);
   [~, curIds] = ismember(curSpecClasses, targetClasses);
    assert(all(curIds));
    
    curSpecAxonIds = cellfun( ...
        @(name) curSpecs.(name).axonIds, ...
        curSpecClasses, 'UniformOutput', false);
    
   [curA, curB] = ndgrid(1:numel(curIds), 1:numel(curIds));
    curSpecMat = zeros(1 + numel(targetClasses));
    curSpecMat(curIds, curIds) = cellfun( ...
        @(idsOne, idsTwo) numel(intersect(idsOne, idsTwo)), ...
        curSpecAxonIds(curA), curSpecAxonIds(curB));
    
    curNonSpecAxonIds = setdiff(curAxonIds, cell2mat(curSpecAxonIds));
    curSpecMat(end, end) = numel(curNonSpecAxonIds);
    
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
    
    plotData(curIdx, :) = curProbs;
end

fig = figure();
fig.Color = 'white';

ax = axes(fig);
plotScaled = plotData ./ sum(plotData, 2);
allBars = bar(ax, plotScaled, 'stacked');

colors = ax.ColorOrder(1:(numel(targetClasses) + 1), :);
mixedColors = (colors(1:(end - 1), :) .^ 2 + colors(2:end, :) .^ 2) / 2;
mixedColors = cat(1, sqrt(mixedColors), zeros(1, 3));

colors = cat(1, transpose(colors), transpose(mixedColors));
colors = transpose(reshape(colors, 3, []));
colors = num2cell(colors(1:(end - 1), :), 2);

[allBars.FaceColor] = deal(colors{:});
[allBars.EdgeColor] = deal('none');
allBars(end).FaceColor = 'none';

xlim(ax, [0.25, numel(plotClasses) + 0.75]);
xticklabels(ax, {axonClasses(plotClasses).tag});
ylabel('Fraction of specific axons');

legends = targetLabels;
legends{end + 1} = 'None';

leg = legend( ...
    allBars(1:2:end), legends, ...
    'Location', 'EastOutside');
leg.Box = 'off';

ax.Box = 'off';
ax.TickDir = 'out';

title( ...
    ax, {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
