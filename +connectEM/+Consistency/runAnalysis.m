% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

thisDir = fileparts(mfilename('fullpath'));
ctrlDir = fullfile(thisDir, 'annotations', 'random-spine-synapses');
twoSpineSynDir = fullfile(thisDir, 'annotations', '2-spine-synapses');

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
ctrlFiles = dir(fullfile(ctrlDir, '*.nml'));
ctrlFiles = fullfile(ctrlDir, {ctrlFiles.name});

ctrlSynT = connectEM.Consistency.loadAnnotations(param, ctrlFiles);
ctrlSynT = vertcat(ctrlSynT{:});

%% Load neurite pairs coupled by two synapses
twoSpineSynFiles = dir(fullfile(twoSpineSynDir, '*.nml'));
twoSpineSynFiles = fullfile(twoSpineSynDir, {twoSpineSynFiles.name});

twoSpineSynT = connectEM.Consistency.loadAnnotations(param, twoSpineSynFiles);
assert(all(cellfun(@(t) size(t, 1), twoSpineSynT) == 2));
twoSpineSynT = vertcat(twoSpineSynT{:});

%% Calculate synapse areas
ctrlSynT.area = ...
    connectEM.Consistency.calcSynapseAreas(param, graph, ctrlSynT);
twoSpineSynT.area = ...
    connectEM.Consistency.calcSynapseAreas(param, graph, twoSpineSynT);

%% Plot synapse area histogram
fig = figure();
fig.Color = 'white';

ax = axes(fig);
hold(ax, 'on');
axis(ax, 'square');
ax.TickDir = 'out';

binEdges = linspace(0, 2, 21);

plotIt = @(values) ...
    histogram( ...
        ax, values, binEdges, ...
        'Normalization', 'probability', ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
plotIt(ctrlSynT.area);
plotIt(twoSpineSynT.area);

xlim(ax, binEdges([1, end]));
xlabel('Axon spine interface area (µm²)');
ylabel('Probability');

legend(ax, ...
    'Random spine synapses', ...
    'Spine synapse pairs', ...
    'Location', 'NorthEast');

title(ax, ...
   {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Plot control consistency
ctrlCVs = combnk(1:numel(ctrlSynT.area), 2);
ctrlCVs = ctrlSynT.area(ctrlCVs);
ctrlCVs = std(ctrlCVs, 0, 2) ./ mean(ctrlCVs, 2);

twoSpineSynCVs = transpose(reshape(twoSpineSynT.area, 2, {}));
twoSpineSynCVs = std(twoSpineSynCVs, 0, 2) ./ mean(twoSpineSynCVs, 2);

fig = figure();
fig.Color = 'white';

ax = axes(fig);
hold(ax, 'on');
axis(ax, 'square');
ax.TickDir = 'out';

binEdges = linspace(0, 1.5, 21);

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

legend(ax, ...
    'Random spine synapses', ...
    'Spine synapse pairs', ...
    'Location', 'NorthEast');

title(ax, ...
   {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
