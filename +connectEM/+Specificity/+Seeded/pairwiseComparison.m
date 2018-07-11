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

curPlot = @(ax, data, varargin) ...
    histogram( ...
        ax, data, ...
        'BinEdges', curBinEdges, ...
        'Normalization', 'probability', ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2, ...
        'FaceAlpha', 1, ...
        varargin{:});

for curAxonClass = reshape(axonClasses, 1, [])
    curAxonIds = curAxonClass.axonIds;
    
    curFig = figure();
    curFig.Color = 'white';
    curFig.Position(3:4) = [1280, 1350];

    for curRow = 1:(curClassCount - 1)
        for curCol = 1:(curClassCount - 1)
           [curSpecsA, curSpecsB] = ...
                pairwiseSeededTargetClassSpec( ...
                    classConn, [curRow, curCol], curAxonIds);
            curPVal = ...
                connectEM.Specificity.kolmogorovSmirnovTest( ...
                    curSpecsA, curSpecsB);

            curAx = curCol + (curRow - 1) * (curClassCount - 1);
            curAx = subplot(curClassCount - 1, curClassCount - 1, curAx);

            hold(curAx, 'on');
            
            curPlot(curAx, curSpecsA(:, curRow));
            curPlot(curAx, curSpecsB(:, curRow));
            
            curSpecs = classConn(curAxonIds, :);
            curSpecs = curSpecs(:, curRow) ./ sum(curSpecs, 2);
            curPlot(curAx, curSpecs, 'EdgeColor', 'black');
            
            xlabel(curAx, { ...
                'Fraction of synapses'; ...
                sprintf('onto %s', targetClasses{curRow})});
    
            annotation( ...
                curFig, ...
                'textbox', curAx.Position, ...
                'String', sprintf('p = %g', curPVal), ...
                'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'top', ...
                'EdgeColor', 'none');
        end
        
        curAxes = flip(curFig.Children(1:(curClassCount - 1)));
        set(curAxes, 'YScale', 'log');
        
        curMaxY = max(arrayfun(@(a) a.YLim(end), curAxes));
        set(curAxes, 'YLim', [0, curMaxY]);
    end

    curAxes = flip(curFig.Children);
    
    for curRow = 1:(curClassCount - 1)
        curAx = curAxes(1 + (curRow - 1) * (curClassCount - 1));
        ylabel(curAx, sprintf('%s-seeded', targetClasses{curRow}));
    end
    
    for curCol = 1:(curClassCount - 1)
        curAx = curAxes(curCol);
        title(curAx, ...
            sprintf('%s-seeded', targetClasses{curCol}), ...
            'FontWeight', 'normal', 'FontSize', 11);
    end
    
    set(curAxes, ...
        'XLim', [0, 1], ...
        'TickDir', 'out', ...
        'PlotBoxAspectRatio', [1, 1, 1], ...
        'DataAspectRatioMode', 'auto');
    
    annotation( ...
        curFig, ...
        'textbox', [0, 0.9, 1, 0.1], ...
        'String', { ...
            info.filename; ...
            info.git_repos{1}.hash; ...
            curAxonClass.title}, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center');
end

%% Utilities
function [specsA, specsB] = ...
        pairwiseSeededTargetClassSpec(classConn, targetIds, axonIds)
    if ~exist('axonIds', 'var') || isempty(axonIds)
        axonIds = reshape(1:size(classConn, 1), [], 1);
    end
    
    specsA = axonTargetClassSpec(classConn, axonIds, targetIds(1));
    specsB = axonTargetClassSpec(classConn, axonIds, targetIds(2));
end

function specs = axonTargetClassSpec(classConn, axonIds, targetId)
    axonIds = axonIds(classConn(axonIds, targetId) > 0);
    
    specs = classConn(axonIds, :);
    specs(:, targetId) = specs(:, targetId) - 1;
    specs = specs ./ sum(specs, 2);
end
