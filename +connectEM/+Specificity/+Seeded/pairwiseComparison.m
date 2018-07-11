% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

minSynPre = 10;

targetClasses = { ...
    'Somata', 'ProximalDendrite', 'ApicalDendrite', ...
    'SmoothDendrite', 'AxonInitialSegment', 'OtherDendrite'};

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, axonClasses] = connectEM.Connectome.load(param, connFile);

[conn, axonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, axonClasses, 'minSynPre', minSynPre);

classConn = ...
    connectEM.Connectome.buildClassConnectome( ...
        conn, 'targetClasses', targetClasses);
    
%% Plot
curClassCount = numel(targetClasses);
curBinEdges = linspace(0, 1, 21);

curPlot = @(ax, data) histogram( ...
    ax, data, ...
    'BinEdges', curBinEdges, ...
    'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);

curFig = figure();
curFig.Color = 'white';

for curRow = 1:(curClassCount - 1)
    for curCol = 1:(curClassCount - 1)
       [curSpecsA, curSpecsB] = ...
            pairwiseSeededTargetClassSpec(classConn, [curRow, curCol]);
        
        curAxIdx = curCol + (curRow - 1) * (curClassCount - 1);
        curAx = subplot(curClassCount - 1, curClassCount - 1, curAxIdx);
        
        hold(curAx, 'on');
        curPlot(curAx, curSpecsA(:, curRow));
        curPlot(curAx, curSpecsB(:, curCol));
    end
end

curAxes = curFig.Children;
set(curAxes, ...
    'XLim', [0, 1], ...
    'TickDir', 'out', ...
    'PlotBoxAspectRatio', [1, 1, 1], ...
    'DataAspectRatioMode', 'auto');

%% Utilities
function [specsA, specsB] = ...
        pairwiseSeededTargetClassSpec(classConn, targetIds)
    axonIdsA = find(classConn(:, targetIds(1)));
    axonIdsB = find(classConn(:, targetIds(2)));
    
    specsA = axonTargetClassSpec(classConn, axonIdsA, targetIds(1));
    specsB = axonTargetClassSpec(classConn, axonIdsB, targetIds(2));
end

function specs = axonTargetClassSpec(classConn, axonIds, targetId)
    specs = classConn(axonIds, :);
    specs(:, targetId) = specs(:, targetId) - 1;
    specs = specs ./ sum(specs, 2);
end
