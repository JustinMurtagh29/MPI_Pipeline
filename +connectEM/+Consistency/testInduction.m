% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

plotConfig = struct;
plotConfig.regionId = 2;
plotConfig.synIds = @(asiT) find( ...
    asiT.type == 'PrimarySpine' ...
  & asiT.axonClass == 'Corticocortical');
plotConfig.title = 'CC primary spine synapses';
plotConfig.tag = 'cc pri sp';

% Plot config
limX = [0, 1.5];
ticksX = linspace(limX(1), limX(2), 4);

limY = [-1.5, 0.5];
ticksY = linspace(limY(1), limY(2), 5);

mapSize = [256, 256];
densityMethod = 'kde2d';

info = Util.runInfo();
Util.showRunInfo(info);

%% Load axon-spine interfaces
[curDir, curAsiFile] = fileparts(connFile);
curAsiFile = fullfile(curDir, sprintf('%s_asiT.mat', curAsiFile));

asiT = load(curAsiFile);
asiT = asiT.asiT;

asiT = asiT(asiT.area > 0, :);
asiT = connectEM.Consistency.Calibration.apply(asiT);

% Finalize plot config
plotConfig.synIds = plotConfig.synIds(asiT);

%% Correlation of synapse size correlation with synapse size
clear cur*;

curPlotConfig = plotConfig;

curSaSdConfig = connectEM.Consistency.buildPairConfigs(asiT, curPlotConfig);
curSaSdConfig = curSaSdConfig(1);
curSaSdConfig.title = curPlotConfig.title;

curCtrlConfig = struct;
curCtrlConfig.synIds = curSaSdConfig.synIdPairs(:);
curCtrlConfig.title = 'SASD';

for i = 1:2
    curKvPairs = { ...
        'xLim', limX, 'yLim', limY, ...
        'mapSize', mapSize, 'method', densityMethod};
    
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
        '%s (n = %d pairs)', curSaSdConfig.title, ...
        size(curSaSdConfig.synIdPairs, 1));
    curCtrlTitle = sprintf( ...
        'vs. random pairs of %s (n = %d)', ...
        curCtrlConfig.title, floor(numel(curCtrlConfig.synIds) / 2));
    
    %% Figure
    curFig = figure();
    curFig.Color = 'white';
    curFig.Position(3:4) = [1060, 970];
    
    curAx = subplot(2, 2, 1);
    image(curAx, ...
        uint8(double(intmax('uint8')) ...
        * curSaSdMap / curMax));
    colormap(curAx, jet(256));
    
    curBar = colorbar('peer', curAx);
    curBar.Ticks = curBar.Limits;
    curBar.TickLabels = {'0', sprintf('%.3g', curMax)};
    
    curAx = subplot(2, 2, 2);
    image(curAx, ...
        uint8(double(intmax('uint8')) ...
        * curCtrlMap / curMax));
    colormap(curAx, jet(256));
    
    curBar = colorbar('peer', curAx);
    curBar.Ticks = curBar.Limits;
    curBar.TickLabels = {'0', sprintf('%.3g', curMax)};
    
    curAx = subplot(2, 2, 3);
    curPValAx = curAx;
    
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
    
    curAx = subplot(2, 2, 4);
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
    
    set( ...
        findobj(curFig.Children, 'Type', 'ColorBar'), ...
        'Location', 'EastOutside', ...
        'TickDirection', 'out', ...
        'Box', 'off');
    
    curAxes = reshape(flip(findobj(curFig, 'Type', 'Axes')), 1, []);
    arrayfun(@(ax) hold(ax, 'on'), curAxes);
    
    curTickIdsX = 1 + floor((mapSize(2) - 1) * ...
        (ticksX - limX(1)) / (limX(2) - limX(1)));
    curTickLabelsX = arrayfun( ...
        @num2str, ticksX, 'UniformOutput', false);
    
    curTickIdsY = 1 + floor((mapSize(1) - 1) * ...
        (ticksY - limY(1)) / (limY(2) - limY(1)));
    curTickLabelsY = arrayfun( ...
        @num2str, ticksY, 'UniformOutput', false);
    
    set(curAxes, ...
        'Box', 'off', ...
        'TickDir', 'out', ...
        'YDir', 'normal', ...
        'YTick', curTickIdsY, ...
        'YTickLabels', curTickLabelsY, ...
        'YLim', [1, mapSize(2)], ...
        'XTick', curTickIdsX, ...
        'XTickLabels', curTickLabelsX, ...
        'XLim', [1, mapSize(1)], ...
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
    
    %% Evaluation    
    curTitle = cell(2, numel(curRegionProps));
    for curRegionId = 1:numel(curRegionProps)
        curMask = curRegionMask == curRegionId;
        curSaSdFrac = sum(curSaSdMap(curMask));
        curDiffFrac = sum(curDiffMap(curMask));
        
        curTitle{1, curRegionId} = ...
            sprintf('%.2f %%', 100 * curSaSdFrac);
        curTitle{2, curRegionId} = ...
            sprintf('%.2f %%', 100 * curDiffFrac);
    end
    
    curTitle = { ...
        'Significance regions'; ...
        strcat(strjoin(curTitle(1, :), ', '), ' of SASD'); ...
        strcat(strjoin(curTitle(2, :), ', '), ' versus control')};
    title(curPValAx, curTitle, 'FontWeight', 'normal', 'FontSize', 10);
    
    %% Build config for next round
    if i > 1; continue; end
    
    % Find synapses in specified region
    curPairT = table;
    curPairT.ids = curSaSdConfig.synIdPairs;
    curPairT.areas = asiT.area(curPairT.ids);
    
    curPairT.x = discretize( ...
        std(curPairT.areas, 0, 2) ./ mean(curPairT.areas, 2), ...
        linspace(limX(1), limX(2), mapSize(2)));
    curPairT.y = discretize( ...
        log10(mean(curPairT.areas, 2)), ...
        linspace(limY(1), limY(2), mapSize(1)));
    
    curPairT.region = sub2ind(mapSize, curPairT.y, curPairT.x);
    curPairT.region = curRegionMask(curPairT.region);
    curMask = curPairT.region == plotConfig.regionId;
    
    % Show size histogram for synapses to be removed
    curBinEdges = linspace(limY(1), limY(2), 21);
    
    curFig = figure();
    curFig.Position(3:4) = [330, 230];
    curAx = axes(curFig); %#ok
    hold(curAx, 'on');
    
    curHist = @(areas) histogram( ...
        curAx, log10(areas), 'BinEdges', curBinEdges);
    
    curHist(curPairT.areas);
    curHist(curPairT.areas(curMask, :));
    
    xlabel(curAx, 'log_{10}(ASI area [µm²])');
    ylabel(curAx, 'Occurences');
    
    connectEM.Figure.config(curFig, info);
    
    % NOTE(amotta): Shuffle synapses outside of specified region. The
    % connections within the specified regions remain as-is.
    rng(0);
    curRandPairIds = curPairT.ids(not(curMask), :);
    curRandPairIds = curRandPairIds(randperm(numel(curRandPairIds)));
    curRandPairIds = reshape(curRandPairIds, [], 2);
    
    curSaSdConfig = struct;
    curSaSdConfig.synIdPairs = [ ...
        curPairT.ids(curMask, :); curRandPairIds];
    curSaSdConfig.title = sprintf( ...
        'Partially shuffled %s', curPlotConfig.title);
end
