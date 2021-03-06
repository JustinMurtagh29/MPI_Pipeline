% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
mix = struct('mean', {}, 'std', {}, 'coeff', {});

i = 1;
% mix(i).mean = -1.00; mix(i).std = 0.2; mix(i).coeff = 0.1; i = i + 1; %#ok
% mix(i).mean = -0.70; mix(i).std = 0.3; mix(i).coeff = 0.8; i = i + 1; %#ok
% mix(i).mean = -0.25; mix(i).std = 0.2; mix(i).coeff = 0.1; i = i + 1; %#ok

% From Stan fit for all corticocortical synapses
mix(i).mean = -0.72; mix(i).std = 0.28; mix(i).coeff = 0.78; i = i + 1; %#ok
mix(i).mean = -0.44; mix(i).std = 0.19; mix(i).coeff = 0.19; i = i + 1; %#ok
mix(i).mean = -0.23; mix(i).std = 0.10; mix(i).coeff = 0.04; i = i + 1; %#ok

clear i;

pairCount = 5290;

mixTitle = arrayfun( ...
    @(m) sprintf( ...
        '%g×N(%g, %g)', m.coeff, m.mean, m.std), ...
    mix, 'UniformOutput', false);
mixTitle = strjoin(mixTitle, ' + ');

info = Util.runInfo();
Util.showRunInfo(info);

%% Analysis of conditional size distribution
curXVec = linspace(-1.5, 0.5, 201);
curTickIds = 1:50:numel(curXVec);

[curCondMap, ~, curSizeDist] = ...
    connectEM.Consistency.Simulation.gmmConditional(mix, 'xVec', curXVec);

curCondMapRel = curCondMap ./ sum(curCondMap, 1);
curCondMapRelDiff = curCondMapRel - (curSizeDist(:) / sum(curSizeDist(:)));

curDensMap = curSizeDist / sum(curSizeDist(:));
curDensMap = curCondMapRel .* reshape(curDensMap, 1, []);

curFig = figure();
curAx = subplot(1, 4, 1);
plot(curAx, curSizeDist, 'Color', 'black');

curAx = subplot(1, 4, 2);
imagesc(curAx, curDensMap);

curAx = subplot(1, 4, 3);
imagesc(curAx, curCondMap);
caxis(curAx, [0, max(curCondMap(:))]);
title(curAx, 'Conditional size distribution');

curAx = subplot(1, 4, 4);
imagesc(curAx, curCondMapRelDiff);
caxis(curAx, [-1, + 1] * max(abs(curCondMapRelDiff(:))));
title(curAx, 'Difference to independent size distribution');

curAllAxes = flip(findobj(curFig, 'Type', 'Axes'));

curHistAx = subplot(1, 4, 1);
curAxes = setdiff(curAllAxes, curHistAx);

for curAx = reshape(curAxes, 1, [])
    hold(curAx, 'on');
    plot(curAx, ...
        [1, size(curCondMap, 2)], ...
        [1, size(curCondMap, 1)], ...
        '--', 'Color', 'white', 'LineWidth', 2);
end

arrayfun(@(ax) axis(ax, 'square'), curAllAxes);
arrayfun(@(ax) colormap(ax, 'jet'), curAllAxes);
set(curAllAxes, 'YDir', 'normal');

set(curAllAxes, ...
    'XLim', [1, size(curCondMap, 2)], ...
    'XTick', curTickIds, 'XTickLabel', arrayfun( ...
    @num2str, curXVec(curTickIds), 'UniformOutput', false));

set(curAxes, ...
    'YLim', [1, size(curCondMap, 1)], ...
    {'YTick'}, get(curAxes, {'XTick'}), ...
    {'YTickLabel'}, get(curAxes, {'XTickLabel'}));

xlabel(curAllAxes(1), 'log10(ASI area [µm²])');
ylabel(curAllAxes(1), 'Probability density');

xlabel(curAllAxes(2), 'log10(ASI area [µm²])');
ylabel(curAllAxes(2), 'log10(ASI area [µm²])');

xlabel(curAllAxes(3), 'log10(ASI area 1 [µm²])');
ylabel(curAllAxes(3), 'log10(ASI area 2 | ASI area 1 [µm²])');

annotation( ...
    curFig, 'textbox', [0, 0.9, 1, 0.1], 'String', mixTitle, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

connectEM.Figure.config(curFig, info);
curFig.Position(3:4) = [1620, 410];

%% Generate synapse pairs
clear cur*;
rng(0);

pairT = table;
[pairT.areas, pairT.classId] = ...
    connectEM.Consistency.Simulation.sampleGmm( ...
        mix, 'numSynapses', 2, 'numConnections', pairCount);

%% Two-dimensional analysis
clear cur*;

% Parameters from connectEM.Connectome.plotSynapseSizeConsistency. The
% bandwidths are the ones returned by `connectEM.Consistency.densityMap`
% for all excitatory connections.
curLimX = [0, 2];
curLimY = [-1.5, 0.5];
curScaleY = 'log10';
curImSize = [256, 256];
curMethod = 'kde2d';

curPvalThreshs = [0.005]; %#ok

curSaSdT = table;
curSaSdT.classId = pairT.classId;
curSaSdT.areas = 10 .^ pairT.areas;
curSaSdT.x = abs( ...
    diff(curSaSdT.areas, 1, 2) ...
    ./ mean(curSaSdT.areas, 2));
curSaSdT.y = log10(mean(curSaSdT.areas, 2));

curSaSdT = curSaSdT( ...
    curLimX(1) <= curSaSdT.x & curSaSdT.x <= curLimX(2) ...
  & curLimY(1) <= curSaSdT.y & curSaSdT.y <= curLimY(2), :);

curSaSdT.mapIds = [ ...
    discretize(curSaSdT.y, linspace( ...
        curLimY(1), curLimY(2), curImSize(1) + 1)), ...
    discretize(curSaSdT.x, linspace( ...
        curLimX(1), curLimX(2), curImSize(2) + 1))];
curSaSdT.mapIdx = sub2ind( ...
    curImSize, curSaSdT.mapIds(:, 1), curSaSdT.mapIds(:, 2));

% Heatmaps
curKvPairs = { ...
    'xLim', curLimX, ...
    'yLim', curLimY, 'yScale', curScaleY, ...
    'mapSize', curImSize, 'method', curMethod};

[curTrueSaSdMap, curBandWidth] = ...
    connectEM.Consistency.densityMap( ...
        curSaSdT.areas, curKvPairs{:});
curTrueRandMaps = ...
    connectEM.Consistency.nullDensityMaps( ...
        curSaSdT.areas, curKvPairs{:}, ...
        'numMaps', 5000, 'bandWidth', curBandWidth);

% Plotting
curMax = max(max(curTrueSaSdMap(:)), max(curTrueSaSdMap(:)));

curTrueRandMap = mean(curTrueRandMaps, 3);
curTrueDiffMap = curTrueSaSdMap - curTrueRandMap;
curMaxDiff = max(abs(curTrueDiffMap(:)));

curPvalMap = 1 - mean(curTrueRandMaps < curTrueSaSdMap, 3);

curPvalImg = -log10(min( ...
    1 - mean(curTrueRandMaps < curTrueSaSdMap, 3), ...
    1 - mean(curTrueRandMaps > curTrueSaSdMap, 3)));

% NOTE(amotta): Detect statistically significant regions. Drop tiny
% regions (with less than 100 pixels or less than 1 % of SASD
% connections), which are most likely caused by outliers.
curRegionMask = curPvalMap < curPvalThreshs(end);
curRegionMask = bwlabel(curRegionMask);

curRegionProps = regionprops( ...
    curRegionMask, {'Area', 'Centroid'}); %#ok
curKeepRegionIds = find([curRegionProps.Area] >= 100);

curKeepRegionIds = curKeepRegionIds(arrayfun( ...
    @(id) sum(curTrueSaSdMap(curRegionMask(:) == id)), ...
    curKeepRegionIds) > 0.01);

curRegionProps = curRegionProps(curKeepRegionIds);
[~, curRegionMask] = ismember(curRegionMask, curKeepRegionIds);
curSaSdT.regionId = curRegionMask(curSaSdT.mapIdx);

%% Figure
curFig = figure();

curAx = subplot(2, 3, 1);
imagesc(curAx, curTrueSaSdMap);
caxis(curAx, [0, curMax]);
colormap(curAx, jet(256));

curAx = subplot(2, 3, 2);
imagesc(curAx, curTrueRandMap);
caxis(curAx, [0, curMax]);
colormap(curAx, jet(256));

curAx = subplot(2, 3, 4);
imagesc(curAx, curPvalImg);
colormap(curAx, jet(256));

hold(curAx, 'on');
for curPvalThresh = curPvalThreshs
    contour(curAx, ...
        curRegionMask & ...
        curPvalMap < curPvalThresh, ...
        true, 'LineColor', 'black');
end

for curRegionId = 1:numel(curRegionProps)
    curPos = curRegionProps(curRegionId).Centroid;

    text(curAx, ...
        curPos(1), curPos(2), num2str(curRegionId), ...
        'Color', 'white', 'FontWeight', 'bold');
end

curTitle = cell(numel(curRegionProps), 1);
for curRegionId = 1:numel(curRegionProps)
    curFracs = nan(numel(curPvalThreshs), 2);
    for curPvalIdx = 1:numel(curPvalThreshs)
        curPvalThresh = curPvalThreshs(curPvalIdx);

        curMask = curRegionMask == curRegionId;
        curMask = curMask & (curPvalMap < curPvalThresh);

        curFracs(curPvalIdx, 1) = sum(curTrueSaSdMap(curMask));
        curFracs(curPvalIdx, 2) = sum(curTrueRandMap(curMask));
    end

    curTitle{curRegionId} = sprintf( ...
        'Region %d: %s', curRegionId, strjoin(arrayfun( ...
            @(v) num2str(v, '%.1f%%'), 100 * curFracs(:, 1), ...
            'UniformOutput', false), ', '));
end

curTitle = [{'Significance regions'}; curTitle];
title(curAx, curTitle, 'FontWeight', 'normal', 'FontSize', 10);

curAx = subplot(2, 3, 5);
imagesc(curAx, curTrueSaSdMap - curTrueRandMap);
caxis(curAx, [-1, +1] * curMaxDiff);
colormap(curAx, jet(256));

curAx = subplot(2, 3, 6);
hold(curAx, 'on');

for curPvalThresh = curPvalThreshs
    contour(curAx, ...
        curRegionMask & ...
        curPvalMap < curPvalThresh, ...
        true, 'LineColor', 'black');
end

for curRegionId = 1:numel(curRegionProps)
    curPos = curRegionProps(curRegionId).Centroid;

    text(curAx, ...
        curPos(1), curPos(2), num2str(curRegionId), ...
        'Color', 'white', 'FontWeight', 'bold');
end

[~, curMinGaussIds] = sort([mix.coeff], 'ascend');
curMinGaussIds = curMinGaussIds(1:(end - 1));

for curIdx = 1:numel(curMinGaussIds)
    curColor = get(groot, 'defaultAxesColorOrder');
    curColor = curColor(curIdx, :);
    
    curGaussId = curMinGaussIds(curIdx);
    curMask = curSaSdT.classId == curGaussId;
    
    scatter(curAx, ...
        curSaSdT.mapIds(curMask, 2), ...
        curSaSdT.mapIds(curMask, 1), ...
        '.', 'MarkerEdgeColor', curColor);
end

curAxes = findobj(curFig, 'Type', 'Axes');

set(curAxes, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'YDir', 'normal', ...
    'YLim', [1, curImSize(2)], ...
    'XLim', [1, curImSize(1)], ...
    'PlotBoxAspectRatio', [1, 1, 1], ...
    'DataAspectRatioMode', 'auto');

curTitle = annotation( ...
    curFig, 'textbox', [0, 0.9, 1, 0.1], 'String', mixTitle, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');
connectEM.Figure.config(curFig, info);
curFig.Position(3:4) = [1350, 850];
