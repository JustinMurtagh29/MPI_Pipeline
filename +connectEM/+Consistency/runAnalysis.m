% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

thisDir = fileparts(mfilename('fullpath'));
ctrlDir = fullfile(thisDir, 'annotations', 'random-spine-synapses');

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
ctrlFiles = dir(fullfile(ctrlDir, '*.nml'));
ctrlFiles = fullfile(ctrlDir, {ctrlFiles.name});

ctrlSynT = connectEM.Consistency.loadAnnotations(param, ctrlFiles);
ctrlSynT = vertcat(ctrlSynT{:});

%% Calculate and plot synapse area
ctrlSynT.area = ...
    connectEM.Consistency.calcSynapseAreas( ...
        param, graph, ctrlSynT);

fig = figure();
fig.Color = 'white';

ax = axes(fig);
hold(ax, 'on');
ax.TickDir = 'out';

binEdges = linspace(0, 1.95, 51);

histogram( ...
    ax, ctrlSynT.area, binEdges, ...
    'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);
xlim(ax, binEdges([1, end]));

xlabel('Axon spine interface area (µm²)');
ylabel('Probability');

%% Plot control consistency
ctrlCVs = combnk(1:numel(ctrlSynT.area), 2);
ctrlCVs = ctrlSynT.area(ctrlCVs);
ctrlCVs = std(ctrlCVs, 0, 2) ./ mean(ctrlCVs, 2);

fig = figure();
fig.Color = 'white';

ax = axes(fig);
hold(ax, 'on');
ax.TickDir = 'out';

binEdges = linspace(0, 1.5, 21);

histogram( ...
    ax, ctrlCVs, binEdges, ...
    'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2)
xlim(ax, binEdges([1, end]));

xlabel('Coefficient of variation');
ylabel('Probability');
