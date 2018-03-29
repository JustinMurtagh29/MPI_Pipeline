% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

thisDir = fileparts(mfilename('fullpath'));
ctrlDir = fullfile(thisDir, 'annotations', 'random-spine-synapses');
twoSpineSynDir = fullfile(thisDir, 'annotations', '2-spine-synapses');
fourSpineSynDir = fullfile(thisDir, 'annotations', '4-spine-synapses');

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

%% Loading augmented graph
graph = Graph.load(rootDir);
graph(~graph.borderIdx, :) = [];

borderAreas = fullfile(rootDir, 'globalBorder.mat');
borderAreas = load(borderAreas, 'borderArea2');
borderAreas = borderAreas.borderArea2;

graph.borderArea = borderAreas(graph.borderIdx);
clear borderAreas;

%% Loading control data
loadClosure = @(n) loadSynapses(param, graph, n);
[ctrlSynT, ctrlSynAreas] = loadClosure(ctrlDir);
[twoSpineSynT, twoSpineSynAreas] = loadClosure(twoSpineSynDir);
[fourSpineSynT, fourSpineSynAreas] = loadClosure(fourSpineSynDir);

%% Plot synapse area histogram
legends = { ...
    'Random synapses', ...
    'Synapse pairs', ...
    'Synapse quadruples'};
binEdges = linspace(0, 1.5, 11);

fig = figure();
fig.Color = 'white';

ax = axes(fig);
hold(ax, 'on');
axis(ax, 'square');
ax.TickDir = 'out';

plotIt = @(values) ...
    histogram( ...
        ax, values, binEdges, ...
        'Normalization', 'probability', ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
plotIt(ctrlSynT.area);
plotIt(twoSpineSynT.area);
plotIt(fourSpineSynT.area);

xlim(ax, binEdges([1, end]));
xlabel(ax, 'Axon spine interface area (µm²)');
ylabel(ax, 'Probability');

legend(ax, ...
    legends, 'Location', 'NorthEast');
title(ax, ...
   {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Box plot of synapse areas
fig = figure();
fig.Color = 'white';

ax = axes(fig);
hold(ax, 'on');
axis(ax, 'square');
ax.TickDir = 'out';

groupedAreas = { ...
    ctrlSynT.area; ...
    twoSpineSynT.area; ...
    fourSpineSynT.area};
groupVar = repelem( ...
    reshape(1:numel(groupedAreas), [], 1), ...
    cellfun(@numel, groupedAreas));
groupedAreas = cell2mat(groupedAreas);

boxplot( ...
    ax, groupedAreas, groupVar, ...
    'Labels', legends);

ylim(ax, binEdges([1, end]));
ylabel(ax, 'Axon spine interface area (µm²)');

title(ax, ...
   {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Plot CVs (observed and expected) for all combinations
plotAreas = {twoSpineSynAreas; fourSpineSynAreas};
plotCount = numel(plotAreas);

fig = figure();
fig.Color = 'white';

for curIdx = 1:plotCount
    curAreas = plotAreas{curIdx};
    curAreas = transpose(curAreas);
    curAreas = sort(curAreas, 1, 'descend');
    
    curCoupling = size(curAreas, 1);
    curPairs = sortrows(combnk(1:curCoupling, 2));
   [~, curSortIds] = sort(diff(curPairs, 1, 2), 'ascend');
    curPairs = curPairs(curSortIds, :);
    
    curCVs = reshape(transpose(curPairs), [], 1);
    curCVs = reshape(curAreas(curCVs, :), 2, size(curPairs, 1), []);
    curCVs = squeeze(std(curCVs, 0, 1) ./ mean(curCVs, 1));
    curCVs = transpose(squeeze(curCVs));
    
    curAx = subplot(plotCount, 1, curIdx);
    curLabels = arrayfun( ...
        @(a, b) sprintf('(%d, %d)', a, b), ...
        curPairs(:, 1), curPairs(:, 2), ...
        'UniformOutput', false);
    boxplot(curAx, curCVs, 'Labels', curLabels);
    
    % This needs to come after `boxplot` for some reason.
    curAx.TickDir = 'out';
    curAx.YLim = [0, sqrt(2)];
end

ylabel(curAx, 'Coefficient of variation');

%% Utilities
function [synT, synAreas] = loadSynapses(param, graph, nmlDir)
    import connectEM.Consistency.loadAnnotations;
    import connectEM.Consistency.calcSynapseAreas;
    
    nmlFiles = dir(fullfile(nmlDir, '*.nml'));
    nmlFiles = fullfile(nmlDir, {nmlFiles.name});
    nmlFiles = reshape(nmlFiles, [], 1);
    
    synT = loadAnnotations(param, nmlFiles);
    
    synCount = cellfun(@(t) size(t, 1), synT);
    synCount = unique(synCount);
    
    if isscalar(synCount)
        groupSize = synCount;
    else
        groupSize = 1;
    end
    
    synT = vertcat(synT{:});
    synT.area = calcSynapseAreas(param, graph, synT);
    synAreas = transpose(reshape(synT.area, groupSize, []));
end