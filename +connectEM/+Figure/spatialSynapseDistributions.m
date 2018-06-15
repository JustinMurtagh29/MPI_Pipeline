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
curDimIds = [3, 1];
curMaxImSize = 300;
curBinSizeUm = 2;
curBandWidth = repelem(0.2, 2);

curHalfBoxUm = param.bbox(:, 2) - mean(param.bbox, 2);
curHalfBoxUm = curHalfBoxUm .* param.raw.voxelSize(:) / 1E3;

curLimits = curBinSizeUm * ceil(curHalfBoxUm / curBinSizeUm);
curLimX = [-1, +1] .* curLimits(curDimIds(1));
curLimY = [-1, +1] .* curLimits(curDimIds(2));

curSupport = curHalfBoxUm(flip(curDimIds));
curSupport = transpose(curSupport(:) .* [-1, +1]);

%%
curImSize = [curLimY(2), curLimX(2)];
curImSize = round(curMaxImSize * curImSize / max(curImSize));

[curImGridY, curImGridX] = ndgrid( ...
    linspace(curLimY(1), curLimY(2), curImSize(1)), ...
    linspace(curLimX(1), curLimX(2), curImSize(2)));
curImGrid = cat(2, curImGridY(:), curImGridX(:));

% Precompute images for faster prototyping
curTypeDensities = cell(size(synTypes));
for curIdx = 1:numel(synTypes)
    curSynT = synT(synT.synType == curIdx, :);
    curDataX = curSynT.pos(:, curDimIds(1));
    curDataY = curSynT.pos(:, curDimIds(2));
    
    curIm = ksdensity( ...
       cat(2, curDataY, curDataX), curImGrid, ...
       'Support', curSupport, 'BandWidth', curBandWidth);
    curIm = reshape(curIm, curImSize);
    
    curTypeDensities{curIdx} = curIm;
end

%% Do the actual plotting
curTypeIdx = 2;

curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [720, 480];

curIm = curTypeDensities{curTypeIdx};
curIm = uint8(double(intmax('uint8')) * curIm / max(curIm(:)));

curAx = axes(curFig);
curIm = image(curAx, curIm);
colormap(curAx, jet(256));
axis(curAx, 'image');

curAx.XTick = [];
curAx.YTick = curAx.YTick([1, end]);
curAx.YTickLabel = {'Pia', 'WM'};
curPos = curAx.Position;


plotHist = @(ax, edges, data) ...
    histogram(ax, ...
        'BinEdges', edges, ...
        'BinCounts', data, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2, ...
        'FaceAlpha', 1);


% Synapse histogram along X axis of plot
curEdges = curLimX(1):curBinSizeUm:curLimX(2);
curPosId = synT.pos(:, curDimIds(1));
curPosId = discretize(curPosId, curEdges);
curSynCounts = accumarray( ...
    cat(2, curPosId, synT.synType), ...
    1, [numel(curEdges) - 1, numel(synTypes)]);

curAx = axes(curFig);
plotHist(curAx, curEdges, curSynCounts(:, curTypeIdx));

curAx.YDir = 'reverse';

curAx.Box = 'off';
curAx.TickDir = 'out';
curAx.XAxisLocation = 'top';
curAx.XLim = curEdges([1, end]);
curAx.XTick = [];

curAx.Position(1) = curPos(1);
curAx.Position(3) = curPos(3);
curAx.Position(2) = 0.02;
curAx.Position(4) = curPos(2) - curAx.Position(2);


% Synapse histogram along Y axis of plot
curEdges = curLimY(1):curBinSizeUm:curLimY(2);
curPosId = synT.pos(:, curDimIds(2));
curPosId = discretize(curPosId, curEdges);
curSynCounts = accumarray( ...
    cat(2, curPosId, synT.synType), ...
    1, [numel(curEdges) - 1, numel(synTypes)]);

curAx = axes(curFig);
curHist = plotHist(curAx, curEdges, curSynCounts(:, curTypeIdx));

curHist.Orientation = 'horizontal';
curAx.YDir = 'reverse';

curAx.Box = 'off';
curAx.TickDir = 'out';
curAx.XAxisLocation = 'top';
curAx.YTick = [];
curAx.YLim = curEdges([1, end]);

curAx.Position(1) = curPos(1) + curPos(3);
curAx.Position(3) = 0.98 - curAx.Position(1);
curAx.Position(2) = curPos(2);


annotation( ...
    curFig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
	'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

%% Plotting histogram
%{
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
%}
