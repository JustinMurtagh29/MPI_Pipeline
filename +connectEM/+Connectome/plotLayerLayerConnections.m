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

%% Plot synapse sizes
clear cur*;

curTcL4Mask =  ...
    conn.axonMeta.axonClass(asiT.preAggloId) == 'Thalamocortical' ...
  & conn.denMeta.targetClass(asiT.postAggloId) == 'ProximalDendrite';
curTcAdMask = ...
    conn.axonMeta.axonClass(asiT.preAggloId) == 'Thalamocortical' ...
  & conn.denMeta.targetClass(asiT.postAggloId) == 'ApicalDendrite';
curL4L4Mask = ...
    l4Conn.denMeta.targetClass(l4AsiT.postAggloId) == 'ProximalDendrite';
curL4AdMask = ...
    l4Conn.denMeta.targetClass(l4AsiT.postAggloId) == 'ApicalDendrite';

curPlotT = table;
curPlotT.name = {...
    'L4 → PD'; 'L4 → AD'; ...
    'TC → PD'; 'TC → AD'};
curPlotT.asiAreas = { ...
    l4AsiT.area(curL4L4Mask); ...
    l4AsiT.area(curL4AdMask); ...
    asiT.area(curTcL4Mask); ...
    asiT.area(curTcAdMask)};

curBinEdges = linspace(-1.5, +0.5, 21);

curFig = figure();
curAx = axes(curFig);


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
    hold(ax, 'on');

    cellfun( ...
        @(a) histogram( ...
            ax, log10(a), binEdges, ...
            'Normalization', 'probability', ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2), ...
        plotT.asiAreas(plotIds));

    curLeg = legend(ax, plotT.name(plotIds));
    set(curLeg, 'Location', 'NorthEast', 'Box', 'off');
end
