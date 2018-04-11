% Plots the path length distribution for 
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
trunkFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_v2.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_spine_attachment.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');

targetClasses = { ...
    'AxonInitialSegment', 'AIS'; ...
    'ApicalDendrite', 'AD'; ...
    'SmoothDendrite', 'SD'};

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

% Spine heads
shAgglos = load(shFile);
shAgglos = shAgglos.shAgglos;

% Connectome (and dendrites after spine attachment)
conn = load(connFile);
dendrites = conn.dendrites;

somaAgglos = conn.denMeta.targetClass == 'Somata';
somaAgglos = conn.dendrites(somaAgglos);

%% Calculate path lengths
Util.log('Calculating path length');

dendMeta = conn.denMeta;
[~, dendMeta.targetClassIdx] = ismember( ...
    conn.denMeta.targetClass, targetClasses);
dendMeta(~dendMeta.targetClassIdx, :) = [];

trunkLensUm = ...
    connectEM.Dendrite.calculatePathLengths( ...
        param, dendrites(dendMeta.id), ...
        trunks, shAgglos, somaAgglos);
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
        'LineWidth', 2);
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
