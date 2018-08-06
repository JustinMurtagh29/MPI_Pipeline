% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

binEdges = linspace(-50, 50, 51);
marginUm = 3;

synTypes = { ...
    'Corticocortical', ...
    'Thalamocortical', ...
    'Inhibitory'};

info = Util.runInfo();
Util.showRunInfo(info);

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

%% Numbers
binUm = 5;

curRangeX = param.bbox(1, :) - mean(param.bbox(1, :), 2);
curRangeX = curRangeX .* param.raw.voxelSize(1) / 1E3;

% Correct for Benedikt's synapse margin
curRangeX = curRangeX + marginUm .* [+1, -1];

curBinEdges = [ ...
    -inf, curRangeX(1), curRangeX(1) + binUm, ...
    curRangeX(2) - binUm, curRangeX(2), +inf];
curSynData = accumarray(horzcat( ...
    discretize(synT.pos(:, 1), curBinEdges), synT.synType), ...
    1, [numel(curBinEdges) - 1, numel(synTypes)]);

curSynData = curSynData([1 + 1, end - 1], :);

curInhCount = curSynData(:, 3);
curTcCount = curSynData(:, 2);

curInhExcRatio = curSynData(:, 3) ./ sum(curSynData(:, 1:3), 2);
curTcCcRatio = curSynData(:, 2) ./ sum(curSynData(:, 1:2), 2);

fprintf('Top %d µm of dataset\n', binUm);
fprintf('* Inh synapses: %d\n', curInhCount(1));
fprintf('* TC synapses: %d\n', curTcCount(1));
fprintf('* Inh / (Inh + Exc): %f\n', curInhExcRatio(1));
fprintf('* TC / (TC + CC): %f\n', curTcCcRatio(1));
fprintf('\n');

fprintf('Bottom %d µm of dataset\n', binUm);
fprintf('* Inh synapses: %d\n', curInhCount(end));
fprintf('* TC synapses: %d\n', curTcCount(end));
fprintf('* Inh / (Inh + Exc): %f\n', curInhExcRatio(end));
fprintf('* TC / (TC + CC): %f\n', curTcCcRatio(end));
fprintf('\n');

%% Scatter plot of synpapse position
clear cur*;
curDimIds = [2, 1];
curMaxImSize = 300;
curBinSizeUm = 2;
curBandWidth = repelem(0.2, 2);
curHistSize = 0.1;

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
for curTypeIdx = 1:numel(synTypes)
    curFig = figure();
    curFig.Color = 'white';
    curFig.Position(3:4) = [480, 360];
    curFigRatio = curFig.Position(4) / curFig.Position(3);
    
    curCorrHistSizes = curFig.Position(3:4) / min(curFig.Position(3:4));
    curCorrHistSizes = curHistSize ./ curCorrHistSizes;

    curIm = curTypeDensities{curTypeIdx};
    curIm = uint8(double(intmax('uint8')) * curIm / max(curIm(:)));
    curImRatio = size(curIm, 1) / size(curIm, 2);

    curImAx = axes(curFig); %#ok
    curIm = image(curImAx, curIm);
    colormap(curImAx, jet(256));
    axis(curImAx, 'image');

    curImAx.XTick = [];
    curImAx.YTick = curImAx.YTick([1, end]);
    curImAx.YTickLabel = {'Pia', 'WM'};
    
    curImAx.Position(3:4) = 0.75 * [1, curImRatio / curFigRatio];
    curImAx.Position(1:2) = (1 - curImAx.Position(3:4)) / 2;

    plotHist = @(ax, edges, data) ...
        histogram(ax, ...
            'BinEdges', edges, ...
            'BinCounts', data, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2, ...
            'FaceAlpha', 1);

    % Synapse histogram along Y axis of plot
    curEdges = curLimY(1):curBinSizeUm:curLimY(2);
    curPosId = synT.pos(:, curDimIds(2));
    curPosId = discretize(curPosId, curEdges);
    curSynCounts = accumarray( ...
        cat(2, curPosId, synT.synType), ...
        1, [numel(curEdges) - 1, numel(synTypes)]);

    curHistAx = axes(curFig); %#ok
    curHist = plotHist(curHistAx, curEdges, curSynCounts(:, curTypeIdx));

    curHist.Orientation = 'horizontal';
    curHistAx.YDir = 'reverse';

    curHistAx.Box = 'off';
    curHistAx.TickDir = 'out';
    curHistAx.XAxisLocation = 'top';
    curHistAx.YTick = [];
    curHistAx.YLim = curEdges([1, end]);
    
    curPos = curImAx.Position;
    curHistAx.Position([2, 4]) = curPos([2, 4]);
    curHistAx.Position(1) = curPos(1) + curPos(3);
    curHistAx.Position(3) = curCorrHistSizes(1);

    annotation( ...
        curFig, ...
        'textbox', [0, 0.9, 1, 0.1], ...
        'String', { ...
            info.filename; ...
            info.git_repos{1}.hash; ...
            synTypes{curTypeIdx}}, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center');
end

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

%% Utilities
