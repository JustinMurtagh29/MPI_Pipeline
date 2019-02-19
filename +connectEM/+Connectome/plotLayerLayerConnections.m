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
asiFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified_asiT.mat');
% NOTE(amotta): This file is identical to 20180726T190355_results.mat with
% the exception of the additional `outputMap.axonData.segIds` field.
outputMapFile = '/tmpscratch/amotta/l4/2018-07-26-tracing-based-output-maps/20190117T143833_results.mat';
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');

% NOTE(amotta): Leave empty to avoid debug skeletons
debugDir = '';

asiRunId = '20190211T111321';

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

% Loading tracings and synapses of L4 cells
outputMap = load(outputMapFile);
connFile = outputMap.info.param.connFile;
axonData = outputMap.axonData;

% Load connectome with synapses and TC axons
[conn, syn] = connectEM.Connectome.load(param, connFile);
conn = connectEM.Connectome.prepareForSpecificityAnalysis(conn);

asiT = load(asiFile);
asiT = asiT.asiT;

asiT = asiT(asiT.area > 0, :);
asiT = connectEM.Consistency.Calibration.apply(asiT);

%% Build synapse and axon-spine interface tables
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.axonClass = conn.axonMeta.axonClass(synT.preAggloId);
synT.targetClass = conn.denMeta.targetClass(synT.postAggloId);

[l4Conn, l4SynT, l4AsiT] = forL4( ...
    param, graph, shAgglos, conn, syn, axonData);

l4SynT.targetClass = l4Conn.denMeta.targetClass(l4SynT.postAggloId);
l4AsiT.targetClass = conn.denMeta.targetClass(l4AsiT.postAggloId);

%% Calculate ASI area
clear cur*;

[curDir, curFile] = fileparts(outputMapFile);
curAsiFile = sprintf('%s__%s_asiT.mat', curFile, asiRunId);
curAsiFile = fullfile(curDir, curAsiFile);

if ~exist(curAsiFile, 'file')
    % Calculate position
    curBorderMeta = fullfile(rootDir, 'globalBorder.mat');
    curBorderMeta = load(curBorderMeta, 'borderSize', 'borderCoM');
    curBorderMeta = structfun(@double, curBorderMeta, 'UniformOutput', false);

    % See also connectEM.Consistency.Script.buildAxonSpineInterfaces
    curWmean = @(w, v) ...
        sum((w / sum(w, 1)) .* v, 1);

    l4AsiT.pos = cellfun( ...
        @(ids) curWmean( ...
            curBorderMeta.borderSize(ids), ...
            curBorderMeta.borderCoM(ids, :)), ...
        l4AsiT.borderIds, 'UniformOutput', false);

    l4AsiT.pos = round(cell2mat(l4AsiT.pos));
    l4AsiT.pos(cellfun(@isempty, l4AsiT.borderIds), :) = nan;
    
    % Calculate area
    curIds = find(not(any(isnan(l4AsiT.pos), 2)));
    curAxonAgglos = {axonData.segIds};

    curAreas = nan(numel(curIds), 1);
    parfor curIdx = 1:numel(curIds)
        curId = curIds(curIdx);
        
        curAreas(curIdx) = ...
            connectEM.Consistency.buildAxonSpineInterfaceAreas( ...
                param, curAxonAgglos, shAgglos, l4AsiT(curId, :));
    end
    
    l4AsiT.area = nan(height(l4AsiT), 1);
    l4AsiT.area(curIds) = curAreas;
    Util.save(curAsiFile, l4AsiT, info);
    Util.protect(curAsiFile);
end

l4AsiT = load(curAsiFile);
l4AsiT = l4AsiT.l4AsiT;

l4AsiT = l4AsiT(l4AsiT.area > 0, :);
l4AsiT = connectEM.Consistency.Calibration.apply(l4AsiT);

%% Generate skeleton for debugging
% In particular, I'm interested in L4 → apical dendrite connections. Are
% these true positives? Are they really that small?
clear cur*;

curPoints = Seg.Global.getSegToPointMap(param);
curMst = @(ids, skel) Skeleton.fromMST( ...
    curPoints(ids, :), param.raw.voxelSize, skel);

if ~isempty(debugDir)
    for curAxonId = 1:numel(axonData)
        curAxon = axonData(curAxonId);
       [~, curNmlName] = fileparts(curAxon.nmlFile);
        
        curSynT = l4SynT(l4SynT.preAggloId == curAxonId, :);
        curNumDigits = ceil(log10(1 + height(curSynT)));
        if isempty(curSynT); continue; end
        
        curSynT.segIds = cellfun( ...
            @(a, b) unique(vertcat(a, b)), ...
            syn.synapses.presynId(curSynT.id), ...
            syn.synapses.postsynId(curSynT.id), ...
            'UniformOutput', false);
        
        curSynT.type(:) = categorical({'Shaft'});
        curSynT.asiArea(:) = nan;
        
        curAsiT = l4AsiT(l4AsiT.preAggloId == curAxonId, :);
       [~, curSynId] = ismember(curAsiT.id, curSynT.id);
       
        curSynT.type(curSynId) = curAsiT.type;
        curSynT.asiArea(curSynId) = curAsiT.area;
        
        curSkel = skeleton(curAxon.nmlFile);
        curSkel = Skeleton.setDescriptionFromRunInfo(curSkel, info);
        
        % Get rid of dendrite
        curSkel = curSkel.deleteTreeWithName('Dendrite');
        
        % Add axon agglomerate
        curSkel = curMst(curAxon.segIds, curSkel);
        curSkel.names{end} = 'Axon Agglomerate';        
        
        for curSynIdx = 1:height(curSynT)
            curSyn = curSynT(curSynIdx, :);
            
            curSkel = curMst(curSyn.segIds{1}, curSkel);
            curSkel.names{end} = sprintf( ...
                '%0*d. Synapse %d. %s onto %s. %f µm²', ...
                curNumDigits, curSynIdx, curSyn.id, ...
                curSyn.type, curSyn.targetClass, curSyn.asiArea);
        end
        
        curOutFile = fullfile(debugDir, strcat(curNmlName, '.nml'));
        curSkel.write(curOutFile);
    end
end

%% Analyse TC → L4 and TC → L5 connections
clear cur*;
curConfigs = struct;
curConfigs(1).axonClass = 'Thalamocortical'; curConfigs(1).tag = 'TC';
curConfigs(2).axonClass = 'Corticocortical'; curConfigs(2).tag = 'CC';

for curConfig = curConfigs
    curAsiT = asiT( ...
        asiT.type == 'PrimarySpine' ...
      & asiT.axonClass == curConfig.axonClass, :);

    tcAsiT = table;
    tcAsiT.name = categories(curAsiT.targetClass);
    tcAsiT.name = strcat({[curConfig.tag, ' → ']}, tcAsiT.name);

    tcAsiT.synCount = accumarray(double(curAsiT.targetClass), 1);
    tcAsiT.synFrac = tcAsiT.synCount ./ sum(tcAsiT.synCount);

    curAreas = accumarray( ...
        double(curAsiT.targetClass), ...
        curAsiT.area, [], @(areas) {areas});
    tcAsiT.meanArea = cellfun(@mean, curAreas);
    tcAsiT.medianArea = cellfun(@median, curAreas);

    tcAsiT = tcAsiT(tcAsiT.synCount > 0, :);
    disp(tcAsiT)
end

%% Analyse L4 → L4 and L4 → L5 connections
clear cur*;
curAsiT = l4AsiT(l4AsiT.type == 'PrimarySpine', :);
    
l4Out = table;
l4Out.name = categories(curAsiT.targetClass);
l4Out.name = strcat({'L4 → '}, l4Out.name);

l4Out.synCount = accumarray(double(l4SynT.targetClass), 1);
l4Out.synFrac = l4Out.synCount ./ sum(l4Out.synCount);

curAreas = accumarray( ...
    double(curAsiT.targetClass), ...
    curAsiT.area, [], @(areas) {areas});
l4Out.meanArea = cellfun(@mean, curAreas);
l4Out.medianArea = cellfun(@median, curAreas);

l4Out = l4Out(l4Out.synCount > 0, :);
disp(l4Out);

%% Prepare data for plotting
binEdges = linspace(-1.5, +0.5, 21);

plotT = table;
plotT.name = {...
    'L4 → PD'; 'L4 → AD'; ...
    'TC → PD'; 'TC → AD'};
plotT.asiAreas = { ...
    l4AsiT.area( ...
        l4AsiT.type == 'PrimarySpine' ...
      & l4AsiT.targetClass == 'ProximalDendrite'); ...
    l4AsiT.area( ...
        l4AsiT.type == 'PrimarySpine' ...
      & l4AsiT.targetClass == 'ApicalDendrite'); ...
    asiT.area( ...
        asiT.type == 'PrimarySpine' ...
      & asiT.axonClass == 'Thalamocortical' ...
      & asiT.targetClass == 'WholeCell'); ...
    asiT.area( ...
        asiT.type == 'PrimarySpine' ...
      & asiT.axonClass == 'Thalamocortical' ...
      & asiT.targetClass == 'ApicalDendrite')};
plotT.medianAsiArea = cellfun(@median, plotT.asiAreas);

plotT.name = cellfun(@(name, areas) ...
    sprintf('%s (n = %d)', name, numel(areas)), ...
    plotT.name, plotT.asiAreas, 'UniformOutput', false);
disp(plotT)

%% Histogram of synapse sizes
clear cur*;
curPlotPair = @(ax, ids) plotPair(binEdges, plotT, ax, ids);

curFig = figure();
curAx = subplot(2, 2, 1); curPlotPair(curAx, [1, 2]);
curAx = subplot(2, 2, 2); curPlotPair(curAx, [1, 3]);
curAx = subplot(2, 2, 3); curPlotPair(curAx, [2, 4]);
curAx = subplot(2, 2, 4); curPlotPair(curAx, [3, 4]);

% Cosmetics
curAxes = findobj(curFig, 'type', 'Axes');
curYLims = cell2mat(get(curAxes, {'YLim'}));
curYLims(:, 2) = max(curYLims(:, 2));

set(curAxes, ...
    {'YLim'}, num2cell(curYLims, 2), ...
    'PlotBoxAspectRatio', [1, 1, 1], ...
    'DataAspectRatioMode', 'auto');

curAx = subplot(2, 2, 3);
xlabel(curAx, 'log_{10}(axon-spine interface area [µm²])');
ylabel(curAx, 'Probability');

curFig.Position(3:4) = [650, 540];
connectEM.Figure.config(curFig, info);

%% Box plot of synapse sizes
clear cur*;
curVals = log10(cell2mat(plotT.asiAreas));
curGroups = repelem(1:height(plotT), cellfun(@numel, plotT.asiAreas));

curFig = figure();
curAx = axes(curFig);

boxplot(curVals, curGroups);

% Cosmetics
curFig.Position(3:4) = [320, 422];
ylabel(curAx, 'log_{10}(axon-spine interface area [µm²])');
xticklabels(curAx, plotT.name);
curAx.XTickLabelRotation = 15;
connectEM.Figure.config(curFig, info);

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
            param, graph, shAgglos, l4Conn, syn, 'addBorderIdsVar', true);
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
