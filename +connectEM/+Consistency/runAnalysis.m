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
    'Random spine synapses', ...
    'Spine synapse pairs', ...
    'Spine synapse quadruples'};
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

%% Plot control consistency
legends = { ...
    'Random spine synapses', ...
    'Spine synapse pairs'};
binEdges = linspace(0, 1.5, 16);

ctrlCVs = combnk(1:numel(ctrlSynT.area), 2);
ctrlCVs = ctrlSynT.area(ctrlCVs);
ctrlCVs = std(ctrlCVs, 0, 2) ./ mean(ctrlCVs, 2);

twoSpineSynCVs = transpose(reshape(twoSpineSynT.area, 2, []));
twoSpineSynCVs = std(twoSpineSynCVs, 0, 2) ./ mean(twoSpineSynCVs, 2);

fig = figure();
fig.Color = 'white';

ax = axes(fig);
hold(ax, 'on');
axis(ax, 'square');
ax.TickDir = 'out';

plotIt = @(data) ...
    histogram( ...
        ax, data, binEdges, ...
        'Normalization', 'probability', ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
plotIt(ctrlCVs);
plotIt(twoSpineSynCVs);

xlim(ax, binEdges([1, end]));
xlabel('Coefficient of variation');
ylabel('Probability');

legend(ax, legends, 'Location', 'NorthEast');

title(ax, ...
   {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Box plot of consistency
fig = figure();
fig.Color = 'white';

ax = axes(fig);
hold(ax, 'on');
axis(ax, 'square');
ax.TickDir = 'out';

groupedAreas = { ...
    ctrlCVs; ...
    twoSpineSynCVs};
groupVar = repelem( ...
    reshape(1:numel(groupedAreas), [], 1), ...
    cellfun(@numel, groupedAreas));
groupedAreas = cell2mat(groupedAreas);

boxplot( ...
    ax, groupedAreas, groupVar, ...
    'Labels', legends);

ylim(ax, binEdges([1, end]));
ylabel(ax, 'Coefficient of variation');

title(ax, ...
   {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

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