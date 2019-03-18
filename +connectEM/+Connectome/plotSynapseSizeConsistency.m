% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
asiRunId = '20190227T082543';

info = Util.runInfo();
Util.showRunInfo(info);

%% Load axon-spine interfaces
[curDir, curAsiFile] = fileparts(connFile);
curAsiFile = sprintf('%s__%s_asiT.mat', curAsiFile, asiRunId);
curAsiFile = fullfile(curDir, curAsiFile);

asiT = load(curAsiFile);
asiT = asiT.asiT;

asiT = asiT(asiT.area > 0, :);
asiT = connectEM.Consistency.Calibration.apply(asiT);

%% Prepare data
clear cur*;

plotConfigs = struct('synIds', {}, 'title', {}, 'tag', {});

axonClasses = { ...
    'Exc', {'Corticocortical', 'Thalamocortical'}; ...
    'CC',  {'Corticocortical'}; ...
    'TC',  {'Thalamocortical'}};
targetClasses = { ...
    'All', categories(asiT.targetClass); ...
    'PD',  'ProximalDendrite'; ...
    'AD',  'ApicalDendrite'; ...
    'OD',  'OtherDendrite'};

for curAxonIdx = 1:size(axonClasses, 1)
    curAxonClass = axonClasses(curAxonIdx, :);
    
    for curTargetIdx = 1:size(targetClasses, 1)
        curTargetClass = targetClasses(curTargetIdx, :);
        
        curSynIds = find( ...
            asiT.type == 'PrimarySpine' ...
          & ismember(asiT.axonClass, curAxonClass{2}) ...
          & ismember(asiT.targetClass, curTargetClass{2}));
        if isempty(curSynIds); continue; end
        
        curTitle = sprintf( ...
            '%s → %s primary spine synapses', ...
            curAxonClass{1}, curTargetClass{1});
        curTag = sprintf('%s %s pri sp', ...
            curAxonClass{1}, curTargetClass{1});
        
        curPlotConfig = plotConfigs([]);
        curPlotConfig(1).synIds = curSynIds;
        curPlotConfig(1).title = curTitle;
        curPlotConfig(1).tag = curTag;
        
        plotConfigs(end + 1) = curPlotConfig; %#ok
    end
end

plotConfigs = reshape( ...
    plotConfigs, ...
    size(targetClasses, 1), ...
    size(axonClasses, 1));

%% Set p-value thresholds for evaluation
% NOTE(amotta): p-value threshold at which large and small low-CV regions
% are not yet merged. (I've only looked at range up to p = 10 %.)
curPvalThreshs = [0.5, 1, 2, 3, 4, 5] / 100;

% Excitatory axons
plotConfigs(1, 1).twoDimPValThreshs = 0.026;
plotConfigs(2, 1).twoDimPValThreshs = inf;
plotConfigs(3, 1).twoDimPValThreshs = 0.015;
plotConfigs(4, 1).twoDimPValThreshs = inf;

% Corticocortical axons
plotConfigs(1, 2).twoDimPValThreshs = 0.052;
plotConfigs(2, 2).twoDimPValThreshs = inf;
plotConfigs(3, 2).twoDimPValThreshs = inf;
plotConfigs(4, 2).twoDimPValThreshs = inf;

% Thalamocortical axons
plotConfigs(1, 3).twoDimPValThreshs = 0.087;
plotConfigs(2, 3).twoDimPValThreshs = inf;
plotConfigs(3, 3).twoDimPValThreshs = inf;
plotConfigs(4, 3).twoDimPValThreshs = inf;

for curIdx = 1:numel(plotConfigs)
    curThreshs = plotConfigs(curIdx).twoDimPValThreshs;
    
    % Only consider p-value thresholds before LTP and LTD regions merger
    curThreshs = curPvalThreshs(curPvalThreshs < curThreshs);
    
    % Only consider lowest and highest p-value threshold
    curThreshs = curThreshs(unique([1, numel(curThreshs)]));
    plotConfigs(curIdx).twoDimPValThreshs = curThreshs;
end

%% Report synapse sizes
clear cur*;
fprintf('Synapse sizes\n');

meanAsis = nan(size(plotConfigs));
for curConfigIdx = 1:numel(plotConfigs)
    curConfig = plotConfigs(curConfigIdx);
    curSynT = asiT(curConfig.synIds, :);
    
    curLog10MeanAsi = log10(mean(curSynT.area));
    curMeanLog10Asi = mean(log10(curSynT.area));
    curStdLog10Asi = std(log10(curSynT.area));
    meanAsis(curConfigIdx) = mean(curSynT.area);
    
    fprintf( ...
        '* log10(Mean %s ASI [µm²]): %f\n', ...
        curConfig.title, curLog10MeanAsi);
    fprintf( ...
        '* log10(%s ASI [µm²]): %f ± %f (mean ± std.)\n', ...
        curConfig.title, curMeanLog10Asi, curStdLog10Asi);
    fprintf('\n');
end

curFig = figure();
curAx = axes(curFig);
imagesc(curAx, transpose(meanAsis));
axis(curAx, 'image');

xticks(curAx, 1:size(targetClasses, 1));
xticklabels(curAx, targetClasses(:, 1));

yticks(curAx, 1:size(axonClasses, 1));
yticklabels(curAx, axonClasses(:, 1));

curCbar = colorbar('peer', curAx);
curCbar.Label.String = 'Average ASI area [µm²]';

curFig.Position(3:4) = [280, 180];
connectEM.Figure.config(curFig, info);

%% Plot distribution of synapse size
connectEM.Consistency.plotSizeHistogram( ...
    info, asiT, plotConfigs(1, :), 'scale', 'log10');
connectEM.Consistency.plotSizeHistogram( ...
    info, asiT, plotConfigs(:, 1), 'scale', 'log10');
connectEM.Consistency.plotSizeHistogram( ...
    info, asiT, plotConfigs(2:end, 2:end), 'scale', 'log10');

%% Report number of occurences for degree of coupling
clear cur*;
curPlotConfig = plotConfigs(1, 1);
curSynT = asiT(curPlotConfig.synIds, :);
[~, ~, curDegreeOfCoupling] = unique(curSynT( ...
    :, {'preAggloId', 'postAggloId'}), 'rows');
curDegreeOfCoupling = accumarray(curDegreeOfCoupling, 1);

docT = table;
[docT.degreeOfCoupling, ~, curDocCount] = unique(curDegreeOfCoupling);

docT.occurences = accumarray(curDocCount, 1);
docT(docT.degreeOfCoupling == 1, :) = [];

fprintf('\nDegree of coupling histogram\n');
disp(docT)

%% Plot histogram of degree of coupling
connectEM.Consistency.plotCouplingHistogram( ...
    info, asiT, plotConfigs(1, 1), 'normalization', 'count');

%% Average ASI area distribution
clear cur*;

curBinEdges = linspace(-1.5, 0.5, 21);

for curIdx = 1%:numel(plotConfigs)
    curPlotConfig = plotConfigs(curIdx);
    
    curSaSdConfig = ...
        connectEM.Consistency.buildPairConfigs(asiT, curPlotConfig);
    curSaSdConfig = curSaSdConfig(1);
    
    curCtrlConfig = struct;
    curCtrlConfig.title = 'SASD';
    curCtrlConfig.synIds = unique(curSaSdConfig.synIdPairs(:));
    
    curCtrlConfig = ...
        connectEM.Consistency.buildPairConfigs(asiT, curCtrlConfig);
    curCtrlConfig = curCtrlConfig(end);
    
    curFig = figure();
    curAx = axes(curFig); %#ok
    hold(curAx, 'on');
    
    histogram(curAx, ...
        log10(mean(asiT.area(curSaSdConfig.synIdPairs), 2)), ...
        'BinEdges', curBinEdges, 'Normalization', 'probability');
    histogram(curAx, ...
        log10(mean(asiT.area(curCtrlConfig.synIdPairs), 2)), ...
        'BinEdges', curBinEdges, 'Normalization', 'probability');
    
    xlabel(curAx, 'log10(Average ASI area [µm²])');
    ylabel(curAx, 'Probability');
    
    curLeg = {curCtrlConfig.title, curSaSdConfig.title};
    curLeg = legend(curAx, curLeg, 'Location', 'SouthOutside');
    
    connectEM.Figure.config(curFig, info);
    curFig.Position(3:4) = [300, 260];
end

%% Synapse area variability
clear cur*;
pValues = nan(4, numel(plotConfigs));

for curIdx = 1:numel(plotConfigs)
    curPlotConfig = plotConfigs(curIdx);
    
    curPairConfigs = ...
        connectEM.Consistency.buildPairConfigs(asiT, curPlotConfig);
    
    %{
    % Uncomment to perform control against random pairs from SASD synapses
    curPairConfigs = Util.modifyStruct( ...
        curPlotConfig, 'synIds', unique(curPairConfigs(1).synIdPairs(:)));
    curPairConfigs = ...
        connectEM.Consistency.buildPairConfigs(asiT, curPairConfigs);
    %}
    
    curFig = ...
        connectEM.Consistency.plotVariabilityHistogram( ...
            info, asiT, curPlotConfig, curPairConfigs(:));
    curFig.Position(3:4) = [370, 540];
   [curLearnedFrac, curUnlearnedFrac, curCvThresh] = ...
        connectEM.Consistency.calculateLearnedFraction( ...
            asiT, curPairConfigs(1), curPairConfigs(end), ...
            'method', 'maxdifference');
    curOpenFrac = 1 - curLearnedFrac - curUnlearnedFrac;

    % Quantitative reporting
    fprintf('**%s**\n', curPlotConfig.title);
    fprintf('Relative size difference between pairs (mean ± std)\n');
    
    for curPairConfig = curPairConfigs
        curRelDiff = asiT.area(curPairConfig.synIdPairs);
        curRelDiff = abs(diff(curRelDiff, 1, 2)) ./ mean(curRelDiff, 2);

        fprintf( ...
            '* %s: %f ± %f\n', ...
            curPairConfig.title, ...
            mean(curRelDiff), std(curRelDiff));
    end
    
    fprintf('\n');
    
    fprintf('Areas under and between curves\n');
    fprintf('* Learned fraction: %.1f %%\n', 100 * curLearnedFrac);
    fprintf('* Possibly learned fraction: %.1f %%\n', 100 * curOpenFrac);
    fprintf('* Unlearned fraction: %.1f %%\n', 100 * curUnlearnedFrac);
    fprintf('* Rel. diff. threshold: %.2f\n', sqrt(2) * curCvThresh);
    fprintf('* CV threshold: %.2f\n', curCvThresh);
    fprintf('\n');

    fprintf('Significance tests\n');
    fprintf('p-values for unexpected synapse size similarity\n');
    
    for curPairIdx = 1:4
        curPairConfig = curPairConfigs(curPairIdx);
        
        curPValue = ...
            connectEM.Consistency.testVariability( ...
                asiT, curPairConfig, curPairConfigs(end));
        pValues(curPairIdx, curIdx) = curPValue;
        
        fprintf('* %s: %g\n', curPairConfig.title, curPValue);
    end

    fprintf('\n');
    fprintf('\n');
end

pValues = reshape(pValues, [4, size(plotConfigs)]);

%% Show p-values
clear cur*;

curTitles = ...
    connectEM.Consistency.buildPairConfigs(asiT, plotConfigs(1));
curTitles = curTitles(1:(end - 1));

curTitles = cellfun( ...
    @(n) n(1:(find(n == '(', 1) - 2)), ...
    {curTitles.title}, 'UniformOutput', false);
assert(isequal(numel(curTitles), size(pValues, 1)));

curFig = figure();

for curIdx = 1:numel(curTitles)
    curPValues = pValues(curIdx, :, :);
    curPValues = shiftdim(curPValues, 1);
    curMap = -log10(curPValues);
    
    subplot(2, 2, curIdx);
    imagesc(transpose(curMap));
    title(curTitles{curIdx});
    
    caxis(feval(@(v) [0, v(2)], caxis()));
    colorbar();
end

curAxes = flip(findobj(curFig, 'type', 'axes'));
xticks(curAxes, 1:size(targetClasses, 1));
xticklabels(curAxes, targetClasses(:, 1));

yticks(curAxes, 1:size(axonClasses, 1));
yticklabels(curAxes, axonClasses(:, 1));

curFig.Position(3:4) = [1050, 700];
connectEM.Figure.config(curFig, info);

%% Correlation of synapse size correlation with synapse size
clear cur*;

% Density difference map
curScaleY = 'linear';
curImSize = [256, 256];
curMethod = 'kde2d';

% Show region in which the same-axon same-dendrite primary spine synapse
% pairs with at least one synapse in the smallest 10th percentile are
% located.
curMinPrctiles = 10;

switch lower(curScaleY)
    case 'linear'
        curLimY = [0, 1];
        curFuncY = @(areas) mean(areas, 2);
        curLabelY = 'Average ASI area [µm²]';
    case 'log10'
        curLimY = [-1.5, 0.5];
        curFuncY = @(areas) log10(mean(areas, 2));
        curLabelY = 'log10(Average ASI area [µm²])';
    otherwise
        error('Invalid Y scale "%s"', curScaleY);
end

curLimX = [0, 2];
curTicksX = linspace(curLimX(1), curLimX(2), 5);
curTicksY = linspace(curLimY(1), curLimY(2), 5);

% NOTE(amotta): The `curMinMap` and `curMaxMap` matrices have the same size
% as the heatmaps of the CV × log10(avg. ASI) space. They contain the ASI
% areas of the smaller and larger synapses, respectively.
curLog10Avg = linspace(curLimY(1), curLimY(2), curImSize(1));
if strcmpi(curScaleY, 'linear'); curLog10Avg = log10(curLog10Avg); end
curCv = linspace(curLimX(1), curLimX(2), curImSize(2)) / sqrt(2);
curCv(curCv >= sqrt(2)) = nan;

curMinMap = (10 .^ curLog10Avg(:)) .* (1 - curCv / sqrt(2));
curMaxMap = (10 .^ curLog10Avg(:)) .* (1 + curCv / sqrt(2));

curFig = figure();
curFigs = curFig([]);
close(curFig);

for curConfig = reshape(plotConfigs, 1, [])
    curSaSdConfig = connectEM.Consistency.buildPairConfigs(asiT, curConfig);
    curSaSdConfig = curSaSdConfig(1);
            
    % Connection types
    curSaSdT = table;
    curSaSdT.synIds = curSaSdConfig.synIdPairs;
    curSaSdT.axonClass = asiT.axonClass(curSaSdT.synIds(:, 1));
    curSaSdT.targetClass = asiT.targetClass(curSaSdT.synIds(:, 1));

    curSaSdT.areas = asiT.area(curSaSdT.synIds);
    curSaSdT.x = abs( ...
        diff(curSaSdT.areas, 1, 2) ...
        ./ mean(curSaSdT.areas, 2));
    curSaSdT.y = curFuncY(curSaSdT.areas);
    
    curSaSdT = curSaSdT( ...
        curLimX(1) <= curSaSdT.x & curSaSdT.x <= curLimX(2) ...
      & curLimY(1) <= curSaSdT.y & curSaSdT.y <= curLimY(2), :);

    curSaSdT.mapIdx = [ ...
        discretize(curSaSdT.y, linspace( ...
            curLimY(1), curLimY(2), curImSize(1) + 1)), ...
        discretize(curSaSdT.x, linspace( ...
            curLimX(1), curLimX(2), curImSize(2) + 1))];
    curSaSdT.mapIdx = sub2ind( ...
        curImSize, curSaSdT.mapIdx(:, 1), curSaSdT.mapIdx(:, 2));
    
    curCtrlConfigs = struct('synIds', {}, 'title', {}, 'tag', {});
    curCtrlConfigs(1).synIds = curSaSdConfig.synIdPairs(:);
    curCtrlConfigs(1).title = 'SASD';
    
    for curCtrlConfig = curCtrlConfigs
        curKvPairs = { ...
            'xLim', curLimX, ...
            'yLim', curLimY, 'yScale', curScaleY, ...
            'mapSize', curImSize, 'method', curMethod};
        
       [curSaSdMap, curBw] = ...
            connectEM.Consistency.densityMap( ...
                asiT.area(curSaSdConfig.synIdPairs), curKvPairs{:});
        
        %{
        % NOTE(amotta): Uncomment to simulate a null model in which random
        % synapse pairs are sampled from within the same connection type
        % (i.e., axon and target class pairs.
       [~, ~, curIds] = unique( ...
           curSaSdT(:, {'axonClass', 'targetClass'}), 'rows');
        curSaSdMap = ...
            connectEM.Consistency.nullDensityMaps( ...
                curSaSdT.areas, 'asiGroups', repmat(curIds, 1, 2), ...
                curKvPairs{:}, 'bandWidth', curBw, 'numMaps', 5000);
        curSaSdMap = mean(curSaSdMap, 3);
        %}
        
        curCtrlMaps = ...
            connectEM.Consistency.nullDensityMaps( ...
                asiT.area(curCtrlConfig.synIds), curKvPairs{:}, ...
                'bandWidth', curBw, 'numMaps', 5000);

        %% Prepare for figure
        curCtrlMap = mean(curCtrlMaps, 3);
        
        curMax = max(max(curSaSdMap(:)), max(curCtrlMap(:)));
        curDiffMap = curSaSdMap - curCtrlMap;
        curMaxDiff = max(abs(curDiffMap(:)));
        
        curPvalMap = 1 - mean(curCtrlMaps < curSaSdMap, 3);
        curPvalThreshs = sort(curConfig.twoDimPValThreshs, 'descend');
        
        curPvalImg = -log10(min( ...
            1 - mean(curCtrlMaps < curSaSdMap, 3), ...
            1 - mean(curCtrlMaps > curSaSdMap, 3)));
        
        % NOTE(amotta): Detect statistically significant regions. Drop tiny
        % regions (with less than 100 pixels or less than 1 % of SASD
        % connections), which are most likely caused by outliers.
        curRegionMask = curPvalMap < curPvalThreshs(1);
        curRegionMask = bwlabel(curRegionMask);
        
        curRegionProps = regionprops( ...
            curRegionMask, {'Area', 'Centroid'}); %#ok
        curKeepRegionIds = find([curRegionProps.Area] >= 100);
        
        curKeepRegionIds = curKeepRegionIds(arrayfun( ...
            @(id) sum(curSaSdMap(curRegionMask(:) == id)), ...
            curKeepRegionIds) > 0.01);
        
        curRegionProps = curRegionProps(curKeepRegionIds);
       [~, curRegionMask] = ismember(curRegionMask, curKeepRegionIds);
        curSaSdT.regionId = curRegionMask(curSaSdT.mapIdx);
        
        curConfigTitle = sprintf( ...
            '%s (n = %d pairs)', curConfig.title, ...
            size(curSaSdConfig.synIdPairs, 1));
        curCtrlTitle = sprintf( ...
            'vs. random pairs of %s (n = %d)', ...
            curCtrlConfig.title, floor(numel(curCtrlConfig.synIds) / 2));
        
        %% Figure
        curFig = figure();
        curFig.Color = 'white';
        curFig.Position(3:4) = [1060, 970];

        curAx = subplot(2, 2, 1);
        imagesc(curAx, curSaSdMap);
        caxis(curAx, [0, curMax]);
        colormap(curAx, jet(256));
        
        curBar = colorbar('peer', curAx);
        curBar.Ticks = curBar.Limits;
        curBar.TickLabels = {'0', sprintf('%.3g', curMax)};
        
        curAx = subplot(2, 2, 2);
        imagesc(curAx, curCtrlMap);
        caxis(curAx, [0, curMax]);
        colormap(curAx, jet(256));

        curBar = colorbar('peer', curAx);
        curBar.Ticks = curBar.Limits;
        curBar.TickLabels = {'0', sprintf('%.3g', curMax)};
        
        curAx = subplot(2, 2, 3);
        curPValAx = curAx;
        
        imagesc(curAx, curPvalImg);
        colormap(curAx, jet(256));
        
        curBar = colorbar('peer', curAx);
        curBar.Ticks = curBar.Limits;
        curBar.TickLabels = arrayfun( ...
            @(val) sprintf('%.3g', val), ...
            curBar.Limits, 'UniformOutput', false);
        curBar.Label.String = '-log10(p-value)';
        
        for curRegionId = 1:numel(curRegionProps)
            curPos = curRegionProps(curRegionId).Centroid;
            
            text(curAx, ...
                curPos(1), curPos(2), num2str(curRegionId), ...
                'Color', 'white', 'FontWeight', 'bold');
        end
        
        curAx = subplot(2, 2, 4);
        hold(curAx, 'on');
        
        imagesc(curAx, curDiffMap);
        caxis(curAx, [-1, +1] * curMaxDiff);
        colormap(curAx, jet(256));
        
        for curPvalThresh = curPvalThreshs
            contour(curAx, ...
                curRegionMask & ...
                curPvalMap < curPvalThresh, ...
                true, 'LineColor', 'black');
        end
        
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
        
        for curMinPrctile = curMinPrctiles
            curAreaThresh = prctile( ...
                curSaSdT.areas(:), curMinPrctile);
            
            for curAx = curAxes
                contour(curAx, ...
                    curMinMap > curAreaThresh, true, ...
                    'LineStyle', '--', 'Color', 'white');
            end
        end

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
        
        arrayfun(@(ax) xlabel(ax, ...
            'Relative ASI area difference'), curAxes);
        arrayfun(@(ax) ylabel(ax, curLabelY), curAxes);

        annotation( ...
            curFig, ...
            'textbox', [0, 0.9, 1, 0.1], ...
            'String', { ...
                info.filename; info.git_repos{1}.hash; ...
                curConfigTitle; curCtrlTitle}, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center');
        
        curFigs(end + 1) = curFig; %#ok
        
        %% Evaluation
        curTitle = cell(numel(curRegionProps), 1);
        for curRegionId = 1:numel(curRegionProps)
            curFracs = nan(numel(curPvalThreshs), 2);
            for curPvalIdx = 1:numel(curPvalThreshs)
                curPvalThresh = curPvalThreshs(curPvalIdx);
                
                curMask = curRegionMask == curRegionId;
                curMask = curMask & (curPvalMap < curPvalThresh);
                
                curFracs(curPvalIdx, 1) = sum(curSaSdMap(curMask));
                curFracs(curPvalIdx, 2) = sum(curCtrlMap(curMask));
            end
            
            curTitle{curRegionId} = sprintf( ...
                'Region %d: %s', curRegionId, strjoin(arrayfun( ...
                    @(v) num2str(v, '%.1f%%'), 100 * curFracs(:, 1), ...
                    'UniformOutput', false), ', '));
        end
        
        curTitle = [{'Significance regions'}; curTitle]; %#ok
        title(curPValAx, curTitle, 'FontWeight', 'normal', 'FontSize', 10);
        
        %{
        %% Evaluate pairs with small and large synapses
        curFig = figure();
       
        subplot(2, 3, 1);
        curAx = curFig.Children(1);
        
        curX = curMinMap(:);
       [curX, curIds] = sort(curX, 'ascend', 'MissingPlacement', 'last');
        
        hold(curAx, 'on');
        plot(gca, log10(curX), cumsum(curCtrlMap(curIds)));
        plot(gca, log10(curX), cumsum(curSaSdMap(curIds)));
        
        xlabel(curAx, 'log10(Max ASI area of small synapse [µm²])');
        xlim(curAx, curLimY);
        
        subplot(2, 3, 1 + 3);
        curAx = curFig.Children(1);
        
        % NOTE(amotta): Reuses `curX` from above
        curBinEdges = [-inf; sort(curSaSdT.areas(:)); +inf];
        curX = discretize(curX, curBinEdges);
        curX = curX / (numel(curBinEdges) - 1);
        
        hold(curAx, 'on');
        plot(gca, curX, cumsum(curCtrlMap(curIds)));
        plot(gca, curX, cumsum(curSaSdMap(curIds)));
        
        ylabel(curAx, 'Fraction of synapse pairs');
        xlabel(curAx, 'Max ASI area of small synapse [quantile]');
        xlim(curAx, [0, 1]);
        
        subplot(2, 3, 2);
        curAx = curFig.Children(1);
        
        curX = curMaxMap(:);
       [curX, curIds] = sort(curX, 'descend', 'MissingPlacement', 'last');
       
        hold(curAx, 'on');
        plot(curAx, log10(curX), cumsum(curCtrlMap(curIds)));
        plot(curAx, log10(curX), cumsum(curSaSdMap(curIds)));
        
        xlabel(curAx, 'log10(Min ASI area of large synapse [µm²])');
        xlim(curAx, curLimY);
        
        subplot(2, 3, 2 + 3);
        curAx = curFig.Children(1);
        
        % NOTE(amotta): Reuses `curX` from above
        curBinEdges = [-inf; sort(curSaSdT.areas(:)); +inf];
        curX = discretize(curX, curBinEdges);
        curX = curX / (numel(curBinEdges) - 1);
        
        hold(curAx, 'on');
        plot(gca, curX, cumsum(curCtrlMap(curIds)));
        plot(gca, curX, cumsum(curSaSdMap(curIds)));
        
        xlabel(curAx, 'Min ASI area of large synapse [quantile]');
        xlim(curAx, [0, 1]);
        
        subplot(2, 3, 3);
        curAx = curFig.Children(1);
        
        curX = linspace(curLimY(1), curLimY(2), curImSize(1));
       [curY, curX] = ndgrid(10 .^ curX, 10 .^ curX);
       
        curMap = arrayfun( ...
            @(x, y) sum(curSaSdMap( ...
                curMinMap(:) <= x ...
              & curMaxMap(:) >= y)), ...
            curX, curY);
        
        imagesc(curAx, curMap);
        caxis(curAx, [0, 1]);
        
        curPos = curAx.Position;
        curCbar = colorbar();
        curCbar.Label.String = 'Fraction of synapse pairs';
        curAx.Position = curPos;
        
        curTickIds = 1 + floor((curImSize(1) - 1) * ...
            (curTicksY - curLimY(1)) / (curLimY(2) - curLimY(1)));
        curTickLabels = arrayfun( ...
            @num2str, curTicksY, 'UniformOutput', false);
        
        xticks(curAx, curTickIds); xticklabels(curAx, curTickLabels);
        yticks(curAx, curTickIds); yticklabels(curAx, curTickLabels);
        
        xlabel(curAx, 'log10(Max ASI area of small synapse [µm²])');
        ylabel(curAx, 'log10(Min ASI area of large synapse [µm²])');
        
        curLines = findobj(curFig, 'type', 'line');
        set(curLines, 'LineWidth', 2);
        
        curAxes = flip(findobj(curFig, 'type', 'axes'));
        axis(curAxes, 'xy');
        axis(curAxes, 'square');
        
        curAx = curAxes(2);
        curPos = curAx.Position;
        curLeg = {'Random pairs', 'Same-axon same-dendrite pairs'};
        curLeg = legend(curAx, curLeg);
        curLeg.Location = 'SouthOutside';
        curAx.Position = curPos;
        
        annotation( ...
            curFig, 'textbox', [0, 0.9, 1, 0.1], ...
            'String', {curConfigTitle; curCtrlTitle});
        
        curFig.Position(3:4) = [1100, 770];
        connectEM.Figure.config(curFig, info);
        
        %% Evaluate sensitivity to p-value threshold
        curCtrlCounts = curCtrlMaps > curSaSdMap;
        curCtrlCounts = sum(curCtrlCounts, 3);
        
        curRegOff = 0;
        curRegMap = zeros(size(curCtrlCounts));
        curRegFracs = zeros(0, size(curCtrlMaps, 3));
        
        for curThresh = 1:size(curCtrlMaps, 3)
            curRegs = curCtrlCounts < curThresh;
            curRegs = regionprops(curRegs, 'PixelIdxList');
            curRegs = {curRegs.PixelIdxList};
            
            % Sort by size (for cosmetic reasons)
           [~, curIds] = sort(cellfun(@numel, curRegs));
            curRegs = curRegs(curIds);
            
            curRegIds = cellfun( ...
                @(ids) setdiff(curRegMap(ids), 0), ...
                curRegs, 'UniformOutput', false);
            
            curRegNewIds = find(cellfun(@isempty, curRegIds));
            curRegNewIds = reshape(curRegNewIds, 1, []);
            
            for curRegIdx = curRegNewIds
                curRegOff = curRegOff + 1;
                curRegIds{curRegIdx} = curRegOff;
                
                curRegMap(curRegs{curRegIdx}) = curRegOff;
                curRegFracs(curRegOff, :) = nan;
            end
            
            for curRegIdx = 1:numel(curRegs)
                curRegFrac = sum(curSaSdMap(curRegs{curRegIdx}));
                curRegFracs(curRegIds{curRegIdx}, curThresh) = curRegFrac;
            end
        end
        
        curFig = figure();
        curAx = axes(curFig); %#ok

        curX = linspace(0, 1, size(curCtrlMaps, 3));
        curY = transpose(flip(curRegFracs, 1));
        plot(curAx, curX, curY, 'LineWidth', 2);
        
        xlim(curAx, [0, 0.1]);
        xlabel(curAx, 'p-value threshold');
        ylabel(curAx, 'Fraction of SASD pairs in region');
        title(curAx, curConfigTitle);
        
        curFig.Position(3:4) = [330, 260];
        connectEM.Figure.config(curFig, info);
        %}
    end
end
