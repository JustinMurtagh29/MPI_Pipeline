% This script simulates the synapse-seeded specificity analyses typically
% done manually by our experimentalists.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
param = struct;
param.saveFolder = '/gaba/u/mberning/results/pipeline/20170217_ROI';

connName = 'connectome_axons_18_a_ax_spine_syn_clust';
synFile = fullfile(param.saveFolder, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');

minSynPre = 10;

info = Util.runInfo();

%% Loading data
syn = load(synFile);
conn = connectEM.Connectome.load(param, connName);

%% Build class connectome
[classConn, targetClasses] = ...
    connectEM.Connectome.buildClassConnectome(conn);

%% Simulate specificity analysis for axons seeded at AD shaft synapse
plotTargetClasses = {'ApicalDendrite', 'Somata'};

% Plot
fig = figure;
fig.Color = 'white';
fig.Position(3:4) = [540, 510];

ax = axes(fig);
ax.TickDir = 'out';
axis(ax, 'square');
hold(ax, 'on');

for curIdx = 1:numel(plotTargetClasses)
    curClassName = plotTargetClasses{curIdx};
    
   	obsConn = forTargetClass( ...
        conn, syn, classConn, targetClasses, curClassName);
    
    % Get rid of of axons that have too few synapses
    obsConn(sum(obsConn, 2) < (minSynPre - 1), :) = [];
    obsConn = obsConn ./ sum(obsConn, 2);
    
    curMu = mean(obsConn, 1);
    curSigma = std(obsConn, 0, 1);
    
    e = errorbar( ...
        1:numel(curMu), curMu, curSigma, ...
        'LineWidth', 1.25, ...
        'Marker', '.', ...
        'MarkerSize', 18);
end

xlim(ax, [0.5, (numel(targetClasses) + 0.5)]);
xticks(ax, 1:numel(targetClasses));

labels = arrayfun( ...
    @char, targetClasses, ...
    'UniformOutput', false);
xticklabels(ax, labels);
xtickangle(ax, 30);

ylabel(ax, 'Fractional synapses');
ylim(ax, [0, 1]);

legends = cellfun(@(n) ...
    sprintf('Seeded at %s shaft synapse', n), ...
    plotTargetClasses, 'UniformOutput', false);
legend(legends, 'Location', 'North');

title(ax, ....
   {info.filename; info.git_repos{1}.hash; ...
    sprintf('All axons with ≥ %d synapses', minSynPre)}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Utilities
function obsConn = forTargetClass( ...
        conn, syn, classConn, targetClasses, className)
    % Build list of possible seed synapses
    seedSynT = connectEM.Connectome.buildSynapseTable(conn, syn);
    seedSynT.ontoTargetClass = conn.denMeta.targetClass(seedSynT.postAggloId);

    % Only look at AD shaft synapses
    seedSynT(seedSynT.isSpine, :) = [];
    seedSynT(seedSynT.ontoTargetClass ~= className, :) = [];

    % Calculate observed synapses
    axonIds = unique(seedSynT.preAggloId);
    obsConn = classConn(axonIds, :);
    
    % Remove seed synapse
    classMask = (targetClasses == className);
    obsConn(:, classMask) = obsConn(:, classMask) - 1;
end