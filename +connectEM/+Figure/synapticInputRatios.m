% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

minSyn = 10;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn] = connectEM.Connectome.load(param, connFile);

%% Prepare data
dendMeta = conn.denMeta;
dendMeta(dendMeta.synCount < minSyn, :) = [];

axonClasses = unique(conn.axonMeta.axonClass);

classConn = connectEM.Connectome.buildClassConnectome( ...
    conn, 'axonClasses', axonClasses, 'targetClasses', []);
classConn = transpose(classConn(:, dendMeta.id));

dendMeta.inhRatio = classConn(:, 3) ./ sum(classConn(:, 1:3), 2);
dendMeta.tcRatio = classConn(:, 2) ./ sum(classConn(:, 1:2), 2);

%% Plotting
binEdges = linspace(0, 1, 51);

targetClasses = unique(dendMeta.targetClass);

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [800, 1100];

for curIdx = 1:numel(targetClasses)
    curTargetClass = targetClasses(curIdx);
    curMeta = dendMeta(dendMeta.targetClass == curTargetClass, :);
    
    
    % Inh / (Inh + Exc) ratio
    ax = subplot(numel(targetClasses), 2, (curIdx - 1) * 2 + 1);
    
    histogram(ax, ...
        curMeta.inhRatio, ...
        'BinEdges', binEdges, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
    ylabel(ax, char(curTargetClass));
    
    
    % TC / (TC + CC) ratio
    ax = subplot(numel(targetClasses), 2, (curIdx - 1) * 2 + 2);
    
    histogram(ax, ...
        curMeta.tcRatio, ...
        'BinEdges', binEdges, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
end

axes = flip(fig.Children);
set(axes, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'XLim', binEdges([1, end]));

xlabel(axes(end - 1), 'Inh / (Inh + Exc)');
xlabel(axes(end - 0), 'TC / (TC + CC)');

annotation(fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');
