% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outputMapFile = '/tmpscratch/amotta/l4/2018-07-26-tracing-based-output-maps/20180726T190355_results.mat';

binSizeUm = 20;

% Set output directory to save PNG and EPS files of plots
outputDir = '';

info = Util.runInfo();
Util.showRunInfo(info);

%% Utility functions
plotHist = @(ax, edges, counts, varargin) ...
    histogram(ax, ...
        'BinEdges', edges, ...
        'BinCounts', counts, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2, ...
        'FaceAlpha', 1, ...
        varargin{:});

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

outputMap = load(outputMapFile);
connFile = outputMap.info.param.connFile;
axonData = outputMap.axonData;

[conn, syn] = connectEM.Connectome.load(param, connFile);

%% Build grand-average axon
grandAvgAxon = struct;
grandAvgAxon.nmlFile = 'grand average';
grandAvgAxon.pathLength = nan;

grandAvgAxon.synapses = cat(1, axonData.synapses);
grandAvgAxon.synapses = sortrows(grandAvgAxon.synapses, 'somaDist');

%% Perform statistical test
clear cur*;
curSomaDist = grandAvgAxon.synapses.somaDist;
curIsSpineSyn = syn.isSpineSyn(grandAvgAxon.synapses.id);

[~, curPVal] = ttest2( ...
    curSomaDist(curIsSpineSyn), ...
    curSomaDist(~curIsSpineSyn)) %#ok

%% Determine bin edges
binEdges = max(grandAvgAxon.synapses.somaDist / 1E3);
binEdges = binSizeUm * ceil(binEdges / binSizeUm);
binEdges = 0:binSizeUm:binEdges;

%% Plot axon
clear cur*;
curPlotData = cat(1, axonData(:), grandAvgAxon);
curPlotData = reshape(curPlotData, 1, []);

curNumDigits = ceil(log10(1 + numel(curPlotData)));

for curIdx = 1:numel(curPlotData)
    curAxonData = curPlotData(curIdx);
    curSynapses = curAxonData.synapses;
    if isempty(curSynapses); continue; end
    
    curSynapses.somaDist = curSynapses.somaDist / 1E3;
    curSynapses.ontoSpine = syn.isSpineSyn(curSynapses.id);
   [~, curAxonName] = fileparts(curAxonData.nmlFile);

    curData = discretize( ...
        curSynapses.somaDist, binEdges);
    curData = accumarray( ...
       [curData, 1 + curSynapses.ontoSpine], ...
        1, [numel(binEdges) - 1, 2]);

    curFig = figure();
    curFig.Color = 'white';
    curFig.Position(3:4) = [390, 190];

    curAx = axes(curFig); %#ok
    curAx.TickDir = 'out';
    hold(curAx, 'on');

    curPlotHist = @(counts, varargin) ...
        plotHist(curAx, binEdges, counts, varargin{:});

    curPlotHist(curData(:, 2), 'EdgeColor', 'magenta');
    curPlotHist(curData(:, 1), 'EdgeColor', 'black');
    
    xlabel(curAx, 'Axonal path length to soma (µm)');
    curAx.XLim = binEdges([1, end]);
    curAx.YLim(1) = 0;

    title(curAx, { ...
        info.filename; info.git_repos{1}.hash; ...
        sprintf('Axon %s', curAxonName)}, ...
        'FontWeight', 'normal', 'FontSize', 10);
    
    if ~isempty(outputDir)
        curFigFileName = strrep(curAxonName, ' ', '-');
        curFigFileName = fullfile(outputDir, sprintf( ...
            '%0*d_%s', curNumDigits, curIdx, curFigFileName));
        
        export_fig('-r172', strcat(curFigFileName, '.png'), curFig);
        export_fig('-r172', strcat(curFigFileName, '.eps'), curFig);
    end
end
