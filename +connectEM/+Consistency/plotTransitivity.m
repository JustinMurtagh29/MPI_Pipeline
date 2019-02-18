% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

info = Util.runInfo();
Util.showRunInfo(info);

%% Load axon-spine interfaces
[curDir, curAsiFile] = fileparts(connFile);
curAsiFile = fullfile(curDir, sprintf('%s_asiT.mat', curAsiFile));

asiT = load(curAsiFile);
asiT = asiT.asiT;

asiT = asiT(asiT.area > 0, :);
asiT = connectEM.Consistency.Calibration.apply(asiT);

%% Prepare data
clear cur*;

plotConfig = struct;
plotConfig(1).synIds = find(asiT.type == 'PrimarySpine');
plotConfig(1).title = 'primary spine synapses';
plotConfig(1).tag = 'pri sp';

curPairs = connectEM.Consistency.buildPairConfigs(asiT, plotConfig);
curPairs = curPairs(1).synIdPairs;

pairT = table;
pairT.axonId = asiT.preAggloId(curPairs(:, 1));
pairT.dendId = asiT.postAggloId(curPairs(:, 1));
pairT.synIds = curPairs;

pairT.areas = asiT.area(pairT.synIds);
pairT.log10avgAsi = log10(mean(pairT.areas, 2));
pairT.cv = std(pairT.areas, 0, 2) ./ mean(pairT.areas, 2);

axonT = table;
[axonT.id, ~, pairT.relAxonId] = unique(pairT.axonId);

axonT.pairIds = accumarray( ...
    pairT.relAxonId, ...
    reshape(1:height(pairT), [], 1), ...
    [], @(ids) {reshape(ids, 1, [])});

axonT = axonT(cellfun(@numel, axonT.pairIds) == 2, :);
axonT.pairIds = cat(1, axonT.pairIds{:});

%% Plot size histogram
clear cur*;

curBinEdges = linspace(-1.5, 0.5, 21);

curFig = figure();
curFig.Position(3:4) = [425, 350];
curAx = axes(curFig);
hold(curAx, 'on');

curHistOne = histogram( ...
    curAx, log10(asiT.area(pairT.synIds(axonT.pairIds, :))), ...
    'Normalization', 'probability', ...
    'BinEdges', curBinEdges);
curHistTwo = histogram( ...
    curAx, log10(asiT.area(pairT.synIds)), ...
    'Normalization', 'probability', ...
    'BinEdges', curBinEdges);

axis(curAx, 'square');
xlabel(curAx, 'log_{10}(ASI area)');
ylabel(curAx, 'Probability');

curLeg = legend(curAx, { ...
        sprintf( ...
            'Pairs of SASD pairs (n = %d synapses)', ...
            numel(curHistOne.Data)), ...
        sprintf( ...
            'All SASD pairs (n = %d synapses)', ...
            numel(curHistTwo.Data))}, ...
    'Location', 'SouthOutside');

connectEM.Figure.config(curFig, info);

%% Plot CV histogram
clear cur*;

curBinEdges = linspace(0, 1.5, 16);

curFig = figure();
curFig.Position(3:4) = [450, 350];
curAx = axes(curFig);
hold(curAx, 'on');

curHistOne = histogram( ...
    curAx, pairT.cv(axonT.pairIds), ...
    'Normalization', 'probability', ...
    'BinEdges', curBinEdges);
curHistTwo = histogram( ...
    curAx, pairT.cv, ...
    'Normalization', 'probability', ...
    'BinEdges', curBinEdges);

axis(curAx, 'square');
xlabel(curAx, 'Coefficient of variation');
ylabel(curAx, 'Probability');

legend(curAx, { ...
        sprintf( ...
            'Pairs of SASD pairs (n = %d SASD pairs)', ...
            numel(curHistOne.Data)), ...
        sprintf( ...
            'All SASD pairs (n = %d SASD pairs)', ...
            numel(curHistTwo.Data))}, ...
    'Location', 'SouthOutside');

connectEM.Figure.config(curFig, info);

%% Check relationship between CVs in pairs of SASD pairs
clear cur*;

curCvs = pairT.cv(axonT.pairIds);
curCvs = sort(curCvs, 2, 'ascend');

rng(0);

curRand = randperm(numel(curCvs));
curRand = reshape(curCvs(curRand), [], 2);
curRand = sort(curRand, 2, 'ascend');

curFit = fitlm(curCvs(:, 1), curCvs(:, 2));
curRandFit = fitlm(curRand(:, 1), curRand(:, 2));

curFig = figure();
curFig.Position(3:4) = [400, 390];
curAx = axes(curFig);

hold(curAx, 'on');
axis(curAx, 'square');

curScatter = scatter(curAx, curCvs(:, 1), curCvs(:, 2), '.');
curLineOne = plot(curAx, [0; 1.5], curFit.predict([0; 1.5]));
curLineTwo = plot(curAx, [0; 1.5], curRandFit.predict([0; 1.5]));
set([curLineOne, curLineTwo], 'LineWidth', 2);

xlim(curAx, [0, 1.5]);
ylim(curAx, [0, 1.5]);

xlabel(curAx, 'CV of low-CV pair');
ylabel(curAx, 'CV of high-CV pair');
curLeg = legend(curAx, { ...
    sprintf( ...
        'Pairs of SASD pairs (n = %d pairs)', ...
        numel(curScatter.XData)), ...
    'Linear fit for pairs of SASD pairs', ...
    'Linear fit for random pairs of SASD pairs'});
connectEM.Figure.config(curFig, info);

%% Correlation of synapse size correlation with synapse size
clear cur*;

% Density difference map
curLimX = [0, 1.5];
curLimY = [-1.5, 0.5];
curImSize = [256, 256];
curMethod = 'kde2d';

curTicksX = linspace(curLimX(1), curLimX(2), 4);
curTicksY = linspace(curLimY(1), curLimY(2), 5);

curSaSdConfig = connectEM.Consistency.buildPairConfigs(asiT, plotConfig);
curSaSdConfig = curSaSdConfig(1);

curCtrlConfig = struct('synIds', {}, 'title', {}, 'tag', {});
curCtrlConfig(1).synIds = curSaSdConfig.synIdPairs(:);
curCtrlConfig(1).title = 'SASD';

curKvPairs = { ...
    'xLim', curLimX, 'yLim', curLimY, ...
    'mapSize', curImSize, 'method', curMethod};

[curSaSdMap, curBw] = ...
    connectEM.Consistency.densityMap( ...
        asiT.area(curSaSdConfig.synIdPairs), curKvPairs{:});
curCtrlMaps = ...
    connectEM.Consistency.nullDensityMaps( ...
        asiT.area(curCtrlConfig.synIds), curKvPairs{:}, ...
        'bandWidth', curBw, 'numMaps', 5000);

%% Prepare for figure
curCtrlMap = mean(curCtrlMaps, 3);

curMax = max(max(curSaSdMap(:)), max(curCtrlMap(:)));
curDiffMap = curSaSdMap - curCtrlMap;
curMaxDiff = max(abs(curDiffMap(:)));

% NOTE(amotta): Division by two for Bonferroni correction.
curPvalMap = min( ...
    1 - mean(curCtrlMaps < curSaSdMap, 3), ...
    1 - mean(curCtrlMaps > curSaSdMap, 3)) / 2;
curPvalMap = -log10(curPvalMap);

% NOTE(amotta): Detect statistically significant regions. Drop tiny
% regions (with less than 100 pixels), which are most likely caused
% by outliers.
curRegionMask = curPvalMap > 2.5;
curRegionMask = bwlabel(curRegionMask);

curRegionProps = regionprops( ...
    curRegionMask, {'Area', 'Centroid'}); %#ok
curKeepRegionIds = find([curRegionProps.Area] >= 100);

curRegionProps = curRegionProps(curKeepRegionIds);
[~, curRegionMask] = ismember(curRegionMask, curKeepRegionIds);

curConfigTitle = sprintf( ...
    '%s (n = %d pairs)', plotConfig.title, ...
    size(curSaSdConfig.synIdPairs, 1));
curCtrlTitle = sprintf( ...
    'vs. random pairs of %s (n = %d)', ...
    curCtrlConfig.title, floor(numel(curCtrlConfig.synIds) / 2));

curPairT = pairT;
curPairT.x = ceil( ...
    (pairT.cv - curLimX(1)) ...
    ./ diff(curLimX) * (curImSize(2) - 1));
curPairT.y = ceil( ...
    (pairT.log10avgAsi - curLimY(1)) ...
    ./ diff(curLimY) * (curImSize(1) - 1));

curAxonT = axonT;
curAxonT.regions = nan(height(curAxonT), 2);
curAxonT.regions(:, 1) = ...
    curRegionMask(sub2ind(curImSize, ...
        curPairT.y(curAxonT.pairIds(:, 1)), ...
        curPairT.x(curAxonT.pairIds(:, 1))));
curAxonT.regions(:, 2) = ...
    curRegionMask(sub2ind(curImSize, ...
        curPairT.y(curAxonT.pairIds(:, 2)), ...
        curPairT.x(curAxonT.pairIds(:, 2))));

curRegionOneMap = find(curAxonT.regions == 1);
curRegionOneMap = curAxonT.pairIds(mod( ...
    curRegionOneMap + height(curAxonT), 2 * height(curAxonT)));
curRegionOneMap = connectEM.Consistency.densityMap( ...
    curPairT.areas(curRegionOneMap, :), curKvPairs{:});
    
curRegionTwoMap = find(curAxonT.regions == 2);
curRegionTwoMap = curAxonT.pairIds(mod( ...
    curRegionTwoMap + height(curAxonT), 2 * height(curAxonT)));
curRegionTwoMap = connectEM.Consistency.densityMap( ...
	curPairT.areas(curRegionTwoMap, :), curKvPairs{:});

curRegionMapsMax = max( ...
    max(curRegionOneMap(:)), ...
    max(curRegionTwoMap(:)));

%% Figure

curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [1060, 970];

curAx = subplot(3, 2, 1);
hold(curAx, 'on');
image(curAx, ...
    uint8(double(intmax('uint8')) ...
  * curSaSdMap / curMax));
colormap(curAx, jet(256));

for curId = 1:height(axonT)
    curPair = axonT.pairIds(curId, :);
    
    plot(curAx, ...
        curPairT.x(curPair), ...
        curPairT.y(curPair), ...
        'Color', [1, 1, 1, 0.25], ...
        'Marker', '.');
end

curBar = colorbar('peer', curAx);
curBar.Ticks = curBar.Limits;
curBar.TickLabels = {'0', sprintf('%.3g', curMax)};

curAx = subplot(3, 2, 2);
image(curAx, ...
    uint8(double(intmax('uint8')) ...
  * curCtrlMap / curMax));
colormap(curAx, jet(256));

curBar = colorbar('peer', curAx);
curBar.Ticks = curBar.Limits;
curBar.TickLabels = {'0', sprintf('%.3g', curMax)};

curAx = subplot(3, 2, 3);
imagesc(curAx, curPvalMap);
colormap(curAx, 'jet');

curBar = colorbar('peer', curAx);
curBar.Ticks = curBar.Limits;
curBar.TickLabels = arrayfun( ...
    @(val) sprintf('%.3g', val), ...
    curBar.Limits, 'UniformOutput', false);
curBar.Label.String = '-log_{10}(p-value)';

for curRegionId = 1:numel(curRegionProps)
    curPos = curRegionProps(curRegionId).Centroid;

    text(curAx, ...
        curPos(1), curPos(2), num2str(curRegionId), ...
        'Color', 'white', 'FontWeight', 'bold');
end

curAx = subplot(3, 2, 4);
hold(curAx, 'on');

image(curAx, ...
    uint8(double(intmax('uint8')) ...
  * (1 + curDiffMap / curMaxDiff) / 2));
contour(curAx, curRegionMask, 'LineColor', 'black');
colormap(curAx, jet(256));

curBar = colorbar('peer', curAx);
curBar.Ticks = [ ...
    curBar.Limits(1), ...
    mean(curBar.Limits), ...
    curBar.Limits(end)];
curBar.TickLabels = { ...
    sprintf('%.3g', -curMaxDiff), '0', ...
    sprintf('%.3g', +curMaxDiff)};

curAx = subplot(3, 2, 5);
image(curAx, ...
    uint8(double(intmax('uint8')) ...
  * curRegionOneMap / curRegionMapsMax));
colormap(curAx, jet(256));

curBar = colorbar('peer', curAx);
curBar.Ticks = curBar.Limits;
curBar.TickLabels = {'0', sprintf('%.3g', curMax)};

curAx = subplot(3, 2, 6);
image(curAx, ...
    uint8(double(intmax('uint8')) ...
  * curRegionTwoMap / curRegionMapsMax));
colormap(curAx, jet(256));

curBar = colorbar('peer', curAx);
curBar.Ticks = curBar.Limits;
curBar.TickLabels = {'0', sprintf('%.3g', curMax)};

set( ...
    findobj(curFig.Children, 'Type', 'ColorBar'), ...
    'Location', 'EastOutside', ...
    'TickDirection', 'out', ...
    'Box', 'off');

curAxes = reshape(flip(findobj(curFig, 'Type', 'Axes')), 1, []);
arrayfun(@(ax) hold(ax, 'on'), curAxes);

curTickIdsX = 1 + floor((curImSize(2) - 1) * ...
    (curTicksX - curLimX(1)) / (curLimX(2) - curLimX(1)));
curTickLabelsX = arrayfun( ...
    @num2str, curTicksX, 'UniformOutput', false);

curTickIdsY = 1 + floor((curImSize(1) - 1) * ...
    (curTicksY - curLimY(1)) / (curLimY(2) - curLimY(1)));
curTickLabelsY = arrayfun( ...
    @num2str, curTicksY, 'UniformOutput', false);

set(curAxes, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'YDir', 'normal', ...
    'YTick', curTickIdsY, ...
    'YTickLabels', curTickLabelsY, ...
    'YLim', [1, curImSize(2)], ...
    'XTick', curTickIdsX, ...
    'XTickLabels', curTickLabelsX, ...
    'XLim', [1, curImSize(1)], ...
    'PlotBoxAspectRatio', [1, 1, 1], ...
    'DataAspectRatioMode', 'auto');

arrayfun(@(ax) ylabel(ax, ...
    'log_{10}(Average ASI area [µm²])'), curAxes);
arrayfun(@(ax) xlabel(ax, ...
    'Coefficient of variation'), curAxes);

annotation( ...
    curFig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'String', { ...
        info.filename; info.git_repos{1}.hash; ...
        curConfigTitle; curCtrlTitle}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');
