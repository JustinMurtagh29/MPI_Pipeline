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

%% Add NML file ID (as proxy of cell ID) to synapse tables
for curId = 1:numel(axonData)
    axonData(curId).synapses.nmlId(:) = curId;
end

%% Ignore interneuron axons for PLASS analysis
% NOTE(amotta): It turns out that of the two interneurons in the dataset
% (whole cells 17 and 22), only one has part of its axon in the EM volume.
% But said axon is mostly myelinated and doesn't make a single synapse.
% Still, let's remove these tracings for completeness.

interNeuronGtNmlFiles = { ...
    '5a796fcd67000090172d94f1', ... % Whole cell 17. Ground truth tracing in +connectEM/+WholeCell/+Data/border-cells_axon-dendrites-split/5a796fcd67000090172d94f1.nml
    '5a796fcd67000060172d94e0'};    % Whole cell 22. Ground truth tracing in +connectEM/+WholeCell/+Data/border-cells_axon-dendrites-split/5a796fcd67000060172d94e0.nml

[~, curNmlFiles] = cellfun( ...
    @fileparts, {axonData.nmlFile}, 'UniformOutput', false);
axonData(ismember(curNmlFiles, interNeuronGtNmlFiles)) = [];
clear cur*;

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

axonCounts = arrayfun(...
    @(e) numel(unique( ...
        grandAvgAxon.synapses.nmlId( ...
            grandAvgAxon.synapses.somaDist >= 1E3 * e))), ...
	binEdges);
axonCounts = axonCounts(1:(end - 1));

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
    curFig.Position(3:4) = [400, 340];

    curAx = subplot(2, 1, 1);
    curAx.TickDir = 'out';
    hold(curAx, 'on');

    curPlotHist = @(curAx, counts, varargin) ...
        plotHist(curAx, binEdges, counts, varargin{:});

    curPlotHist(curAx, curData(:, 2), 'EdgeColor', 'magenta');
    curPlotHist(curAx, curData(:, 1), 'EdgeColor', 'black');
    
    xlabel(curAx, 'Axonal path length to soma (µm)');
    curAx.XLim = binEdges([1, end]);
    curAx.YLim(1) = 0;

    title(curAx, { ...
        info.filename; info.git_repos{1}.hash; ...
        sprintf('Axon %s', curAxonName)}, ...
        'FontWeight', 'normal', 'FontSize', 10);
    
    curAx = subplot(2, 1, 2);
    curPlotHist(curAx, axonCounts, 'EdgeColor', 'black');
    
    curAx.YDir = 'reverse';
    curAx.XAxisLocation = 'top';
    
    curAx.XLim = binEdges([1, end]);
    curAx.TickDir = 'out';
    curAx.Box = 'off';
    
    if ~isempty(outputDir)
        curFigFileName = strrep(curAxonName, ' ', '-');
        curFigFileName = fullfile(outputDir, sprintf( ...
            '%0*d_%s', curNumDigits, curIdx, curFigFileName));
        
        export_fig('-r172', strcat(curFigFileName, '.png'), curFig);
        export_fig('-r172', strcat(curFigFileName, '.eps'), curFig);
    end
end
