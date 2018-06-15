% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

binEdges = linspace(-50, 50, 51);

synTypes = { ...
    'Corticocortical', ...
    'Thalamocortical', ...
    'Inhibitory'};

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn] = connectEM.Connectome.load(param, connFile);

%% Preparing data
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.synType = conn.axonMeta.axonClass(synT.preAggloId);

[~, synT.synType] = ismember(synT.synType, synTypes);
synT(~synT.synType, :) = [];

synPos = connectEM.Synapse.calculatePositions(param, syn);
synT.pos = synPos(synT.id, :);
clear synPos;

% Physical units, relative to center
synT.pos = synT.pos - mean(param.bbox, 2)';
synT.pos = synT.pos .* param.raw.voxelSize ./ 1E3;

%% Scatter plot of synpapse position
clear cur*;
curDimIds = [2, 1];
curMaxImSize = 300;
curBandWidth = repelem(0.2, 2);

curLimX = round(( ...
    param.bbox(curDimIds(1), :) ...
    - mean(param.bbox(curDimIds(1), :))) ...
    * param.raw.voxelSize(curDimIds(1)) / 1E3);
curLimY = round(( ...
    param.bbox(curDimIds(2), :) ...
    - mean(param.bbox(curDimIds(2), :))) ...
    * param.raw.voxelSize(curDimIds(2)) / 1E3);

curImSize = [curLimY(2), curLimX(2)];
curImSize = round(curMaxImSize * curImSize / max(curImSize));

[curImGridY, curImGridX] = ndgrid( ...
    linspace(curLimY(1), curLimY(2), curImSize(1)), ...
    linspace(curLimX(1), curLimX(2), curImSize(2)));
curImGrid = cat(2, curImGridY(:), curImGridX(:));

curSupport = cat(2, curLimY(:), curLimX(:));

curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [1200, 325];

for curIdx = 1:numel(synTypes)
    curSynT = synT(synT.synType == curIdx, :);
    curDataX = curSynT.pos(:, curDimIds(1));
    curDataY = curSynT.pos(:, curDimIds(2));
    
    curIm = ksdensity( ...
       cat(2, curDataY, curDataX), curImGrid, ...
       'Support', curSupport, 'BandWidth', curBandWidth);
    curIm = uint8(double(intmax('uint8')) * curIm / max(curIm));
    curIm = reshape(curIm, curImSize);
    
    curAx = subplot(1, numel(synTypes), curIdx);
    imshow(curIm, jet(double(intmax('uint8'))));
    
    title(curAx, ...
        synTypes{curIdx}, ...
        'FontWeight', 'normal', ...
        'FontSize', 10);
end

curAxes = curFig.Children;

set(curAxes, ...
    'Visible', 'on', ...
    'Layer', 'top', ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'XTick', [1, curImSize(2)], ...
    'XTickLabel', {}, ...
    'YTick', [1, curImSize(1)], ...
    'YTickLabel', {});

yticklabels(curAxes(end), {'Pia', 'WM'});

annotation( ...
    curFig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
	'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

%% Plotting data
clear cur*;
curDimLabels = {'X', 'Y', 'Z'};

curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [1500, 675];

plotHist = @(ax, data) ...
    histogram(ax, ...
        'BinEdges', binEdges, ...
        'BinCounts', data, ...
        'Orientation', 'horizontal', ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2, ...
        'FaceAlpha', 1);
    
for curDim = 1:3
    curPosId = synT.pos(:, curDim);
    curPosId = discretize(curPosId, binEdges);
    
    curSynCounts = accumarray( ...
        cat(2, curPosId, synT.synType), ...
        1, [numel(binEdges) - 1, numel(synTypes)]);
    curSynRatios = [ ...
        curSynCounts(:, 3) ./ sum(curSynCounts, 2), ...
        curSynCounts(:, 2) ./ sum(curSynCounts(:, 1:2), 2)];
    
    % Synapse probabilities
    curAx = subplot(2, 3, curDim);
    curAx.YAxis.Direction = 'reverse';
    hold(curAx, 'on');
    
    for curTypeId = 1:size(curSynCounts, 2)
        plotHist(curAx, curSynCounts(:, curTypeId));
    end
    
    % Synapse ratios
    curAx = subplot(2, 3, curDim + 3);
    curAx.YAxis.Direction = 'reverse';
    hold(curAx, 'on');
    
    for curRatioId = 1:size(curSynRatios, 2)
        curRatios = curSynRatios(:, curRatioId);
        curRatios(isnan(curRatios)) = 0;
        plotHist(curAx, curRatios);
    end
    
    ylabel(curAx, sprintf('%s (µm)', curDimLabels{curDim}));
end

curAxes = flip(curFig.Children);

% Make square
[curAxes.PlotBoxAspectRatio] = deal([1, 1, 1]);
[curAxes.DataAspectRatioMode] = deal('auto');

% X labels
xlabel(curAxes(1), 'Probability');
xlabel(curAxes(2), 'Ratio');

set(curAxes(1:2:end), 'XLim', [0, max(cat(2, curAxes(1:2:end).XLim))]);
set(curAxes(2:2:end), 'XLim', [0, 0.4]);

% Legends
curLeg = addLegend(curAxes(end - 1), synTypes, 'Location', 'EastOutside');
curLeg.Box = 'off';

curLeg = addLegend(curAxes(end), ...
   {'Inh / (Inh + Exc)', 'TC / (TC + CC)'}, ...
    'Location', 'EastOutside');
curLeg.Box = 'off';

annotation( ...
    curFig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
	'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

%% Utilities
function leg = addLegend(ax, legends, varargin)
    axPos = ax.Position;
    leg = legend(ax, legends, varargin{:});
    ax.Position = axPos;
end
