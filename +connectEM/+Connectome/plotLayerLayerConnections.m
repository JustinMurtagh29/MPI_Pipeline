% This script analyses {thalamus, L4} → {L4, L5} connections. The proximal
% and apical dendrite target classes serve as proxies for L4 and L5
% dendrites, respectively. The axons we call "thalamocortical" most likely
% originate from VPM (see Bopp et al. 2017 J Neurosci). For L4 axons we use
% the soma-based axon tracings.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
% NOTE(amotta): This file is identical to 20180726T190355_results.mat with
% the exception of the additional `outputMap.axonData.segIds` field.
outputMapFile = '/tmpscratch/amotta/l4/2018-07-26-tracing-based-output-maps/20190117T143833_results.mat';
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);

shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

graph = Graph.load(rootDir);
graph = graph(graph.borderIdx ~= 0, :);

borderAreas = fullfile(rootDir, 'globalBorder.mat');
borderAreas = load(borderAreas, 'borderArea2');
borderAreas = borderAreas.borderArea2;

graph.borderArea = borderAreas(graph.borderIdx);
clear borderAreas;

% Loading tracings and synapses of L4 cells
outputMap = load(outputMapFile);
connFile = outputMap.info.param.connFile;
axonData = outputMap.axonData;

% Load connectome with synapses and TC axons
[conn, syn] = connectEM.Connectome.load(param, connFile);
conn = connectEM.Connectome.prepareForSpecificityAnalysis(conn);

%%
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.axonClass = conn.axonMeta.axonClass(synT.preAggloId);
synT.targetClass = conn.denMeta.targetClass(synT.postAggloId);

asiT = ...
    connectEM.Connectome.buildAxonSpineInterfaces( ...
        param, graph, shAgglos, conn, syn);

[l4Conn, l4SynT, l4AsiT] = ...
    forL4(param, graph, shAgglos, conn, syn, axonData);
l4SynT.targetClass = l4Conn.denMeta.targetClass(l4SynT.postAggloId);
l4AsiT = l4AsiT(l4AsiT.type == 'PrimarySpine', :);

%% Analyse TC → L4 and TC → L5 connections
tcOut = table;
tcOut.name = categories(conn.denMeta.targetClass);
tcOut.name = strcat({'TC → '}, tcOut.name);

curMask = synT.axonClass == 'Thalamocortical';
tcOut.synCount = accumarray(double(synT.targetClass(curMask)), 1);
tcOut.synFrac = tcOut.synCount ./ sum(tcOut.synCount);
disp(tcOut)

%% Analyse L4 → L4 and L4 → L5 connections
l4Out = table;
l4Out.name = categories(l4Conn.denMeta.targetClass);
l4Out.name = strcat({'L4 → '}, l4Out.name);

l4Out.synCount = accumarray(double(l4SynT.targetClass), 1);
l4Out.synFrac = l4Out.synCount ./ sum(l4Out.synCount);
disp(l4Out);

%% Prepare data for plotting
binEdges = linspace(-1.5, +0.5, 21);

plotT = table;
plotT.asiAreas = { ...
    l4AsiT.area(l4Conn.denMeta.targetClass( ...
        l4AsiT.postAggloId) == 'ProximalDendrite'); ...
    l4AsiT.area(l4Conn.denMeta.targetClass( ...
        l4AsiT.postAggloId) == 'ApicalDendrite'); ...
    asiT.area( ...
        conn.axonMeta.axonClass(asiT.preAggloId) == 'Thalamocortical' ...
      & conn.denMeta.targetClass(asiT.postAggloId) == 'ProximalDendrite'); ...
    asiT.area( ...
        conn.axonMeta.axonClass(asiT.preAggloId) == 'Thalamocortical' ...
      & conn.denMeta.targetClass(asiT.postAggloId) == 'ApicalDendrite')};
  
plotT.name = {...
    'L4 → PD'; 'L4 → AD'; ...
    'TC → PD'; 'TC → AD'};
plotT.name = cellfun(@(name, areas) ...
    sprintf('%s (n = %d)', name, numel(areas)), ...
    plotT.name, plotT.asiAreas, 'UniformOutput', false);

%% Histogram of synapse sizes
clear cur*;
curPlotPair = @(ax, ids) plotPair(binEdges, plotT, ax, ids);

curFig = figure();
curAx = subplot(2, 2, 1); curPlotPair(curAx, [1, 2]);
curAx = subplot(2, 2, 2); curPlotPair(curAx, [1, 3]);
curAx = subplot(2, 2, 3); curPlotPair(curAx, [2, 4]);
curAx = subplot(2, 2, 4); curPlotPair(curAx, [3, 4]);

% Cosmetics
curFig.Color = 'white';
curAxes = findobj(curFig, 'type', 'Axes');
curYLims = cell2mat(get(curAxes, {'YLim'}));
curYLims(:, 2) = max(curYLims(:, 2));
set(curAxes, {'YLim'}, num2cell(curYLims, 2));
set(curAxes, 'TickDir', 'out', 'XLim', binEdges([1, end]));

curAx = subplot(2, 2, 3);
xlabel(curAx, 'log_{10}(axon-spine interface area [µm²])');
ylabel(curAx, 'Probability');

annotation(curFig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'String', {info.filename; info.git_repos{1}.hash});

%% Box plot of synapse sizes
clear cur*;
curVals = log10(cell2mat(plotT.asiAreas));
curGroups = repelem(1:height(plotT), cellfun(@numel, plotT.asiAreas));

curFig = figure();
curAx = axes(curFig);

boxplot(curVals, curGroups);

% Cosmetics
curFig.Color = 'white';
curFig.Position(3:4) = [320, 422];
curAx.TickDir = 'out';
curAx.Box = 'off';

ylim(curAx, binEdges([1, end]));
ylabel(curAx, 'log_{10}(axon-spine interface area [µm²])');
xticklabels(curAx, plotT.name);
curAx.XTickLabelRotation = 15;

title(curAx, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Utilities
function [l4Conn, l4SynT, l4AsiT] = forL4( ...
        param, graph, shAgglos, conn, syn, axonData)
    % NOTE(amotta): Because the axons used for this analysis are different
    % from the ones in the dense reconstruction, we have to build a fake
    % connectome.
    clear cur*;

    axonAgglos = reshape({axonData.segIds}, [], 1);
    axonAgglos = Agglo.removeSegmentOverlap(axonAgglos);

    maxSegId = Seg.Global.getMaxSegId(param);
    axonLUT = Agglo.buildLUT(maxSegId, axonAgglos);
    dendLUT = Agglo.buildLUT(maxSegId, conn.dendrites);

    curSyn = syn.synapses;
    curSyn.id = reshape(1:height(curSyn), [], 1);

    curSyn.preAggloId = cellfun( ...
        @(ids) setdiff(axonLUT(ids), 0), ...
        curSyn.presynId, 'UniformOutput', false);
    curSyn = curSyn(cellfun( ...
        @isscalar, curSyn.preAggloId), :);

    curSyn.postAggloId = cellfun( ...
        @(ids) setdiff(dendLUT(ids), 0), ...
        curSyn.postsynId, 'UniformOutput', false);
    curSyn = curSyn(cellfun( ...
        @isscalar, curSyn.postAggloId), :);

    curSyn = [ ...
        cell2mat(curSyn.preAggloId), ...
        cell2mat(curSyn.postAggloId), ...
        curSyn.id];

    curConnectome = table;
   [curConnectome.edges, ~, curIds] = ...
        unique(curSyn(:, 1:2), 'rows');
    curConnectome.synIdx = accumarray( ...
        curIds, curSyn(:, end), [], @(ids) {ids});

    l4Conn = struct;
    l4Conn.axons = axonAgglos;
    l4Conn.denMeta = conn.denMeta;
    l4Conn.dendrites = conn.dendrites;
    l4Conn.connectome = curConnectome;
    clearvars cur* -except curConn;
    
    l4SynT = ...
        connectEM.Connectome.buildSynapseTable(l4Conn, syn);
    l4AsiT = ...
        connectEM.Connectome.buildAxonSpineInterfaces( ...
            param, graph, shAgglos, l4Conn, syn);
end

function plotPair(binEdges, plotT, ax, plotIds)
    colors = get(groot, 'defaultAxesColorOrder');
    colors = num2cell(colors, 2);
    
    hold(ax, 'on');

    cellfun( ...
        @(areas, color) histogram( ...
            ax, log10(areas), binEdges, ...
            'Normalization', 'probability', ...
            'DisplayStyle', 'stairs', ...
            'EdgeColor', color, ...
            'LineWidth', 2, ...
            'FaceAlpha', 1), ...
        plotT.asiAreas(plotIds), colors(plotIds));

    curLeg = legend(ax, plotT.name(plotIds));
    set(curLeg, 'Location', 'NorthEast', 'Box', 'off');
end
