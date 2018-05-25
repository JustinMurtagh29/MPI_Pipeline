% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

minSyn = 50;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = connectEM.Connectome.load(param, connFile);

%% Calculate input ratios
somaMeta = conn.denMeta(conn.denMeta.targetClass == 'Somata', :);
somaMeta(somaMeta.synCount < minSyn, :) = [];

axonClasses = unique(conn.axonMeta.axonClass);

classConn = ...
    connectEM.Connectome.buildClassConnectome( ...
        conn, 'axonClasses', axonClasses, 'targetClasses', []);
classConn = transpose(classConn(:, somaMeta.id));

inhFrac = ...
    classConn(:, axonClasses == 'Inhibitory') ...
 ./ sum(classConn(:, axonClasses ~= 'Other'), 2);

tcExcFrac = sum(classConn(:, ismember( ...
    axonClasses, {'Thalamocortical', 'Corticocortical'})), 2);
tcExcFrac = classConn(:, axonClasses == 'Thalamocortical') ./ tcExcFrac;

%% Plot input ratios
binEdges = linspace(0, 1, 21);
plotHist = @(ax, d) histogram( ...
    ax, d, 'DisplayStyle', 'stairs', ...
    'LineWidth', 2, 'BinEdges', binEdges);

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [725, 380];

ax = subplot(1, 2, 1);
axis(ax, 'square');
hold(ax, 'on');
plotHist(ax, inhFrac(~somaMeta.isInterneuron));
plotHist(ax, inhFrac( somaMeta.isInterneuron));

xlabel(ax, 'I / (I + E)');
ylabel(ax, 'Somata');

ax = subplot(1, 2, 2);
axis(ax, 'square');
hold(ax, 'on');
plotHist(ax, tcExcFrac(~somaMeta.isInterneuron));
plotHist(ax, tcExcFrac( somaMeta.isInterneuron));

xlabel(ax, 'TC / (TC + CC)');

maxY = max(arrayfun(@(a) a.YLim(end), fig.Children));
[fig.Children.YLim] = deal([0, maxY]);

leg = legend( ...
    ax, {'Excitatory cells', 'Interneurons'}, ...
    'Location', 'NorthEast');
leg.Box = 'off';

title = { ...
    info.filename; info.git_repos{1}.hash; sprintf( ...
    'Somata with ≥ %d synapses (n = %d)', minSyn, height(somaMeta))};
annotation( ...
    fig, 'textbox', [0, 0.9, 1, 0.1], 'String', title, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');
