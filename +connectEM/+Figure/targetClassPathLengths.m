% Plots the path length distribution for 
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
trunkFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

minSynPost = 10;

targetClasses = { ...
    'ApicalDendrite', 'AD'; ...
    'SmoothDendrite', 'SD'; ...
    'AxonInitialSegment', 'AIS'};

targetLabels = targetClasses(:, 2);
targetClasses = targetClasses(:, 1);

info = Util.runInfo();

%% Loading data
Util.log('Loading data');
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);

% Dendrite prior to spine attachment
trunks = load(trunkFile);
trunks = trunks.dendrites(trunks.indBigDends);
trunks = Agglo.fromSuperAgglo(trunks);

% Connectome (and dendrites after spine attachment)
conn = load(connFile);
dendrites = conn.dendrites;

%% Calculate path lengths
Util.log('Calculating path length');

dendMeta = conn.denMeta;
[~, dendMeta.targetClassIdx] = ismember( ...
    conn.denMeta.targetClass, targetClasses);
dendMeta(~dendMeta.targetClassIdx, :) = [];

% Ignore targets with less than 10 synapses
dendMeta(dendMeta.synCount < minSynPost, :) = [];

trunkLensUm = ...
    connectEM.Dendrite.calculatePathLengths( ...
        param, dendrites(dendMeta.id), trunks);
trunkLensUm = trunkLensUm / 1E3;

%% Plot histogram
binEdges = linspace(0, 400, 41);

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [400, 240];

ax = axes(fig);
hold(ax, 'on');

for curIdx = 1:numel(targetClasses)
    curLensUm = dendMeta.targetClassIdx == curIdx;
    curLensUm = trunkLensUm(curLensUm);
    
    histogram( ...
        ax, curLensUm, ...
        'BinEdges', binEdges, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2, ...
        'FaceAlpha', 1);
end

ax.TickDir = 'out';
xlim(binEdges([1, end]));
xlabel('Length (µm)');
ylabel('Neurites');

legend( ...
    ax, targetLabels, ...
    'Location', 'NorthEast', ...
    'Box', 'off');

title(ax, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
