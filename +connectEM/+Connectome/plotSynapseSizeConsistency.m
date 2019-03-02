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

%% Fit mixture of Gaussians to size distribution
clear cur*;
rng(0);

curN = 2;
curConfig = plotConfigs(2, 1);
curBinEdges = linspace(-2, 0.5, 51);

curLog10Areas = log10(asiT.area(curConfig.synIds));
curGmm = fitgmdist(curLog10Areas, curN);

curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');

% GMM components
curLegs = cell(curN, 1);
for curId = 1:curN
    curMean = curGmm.mu(curId);
    curStd = sqrt(curGmm.Sigma(curId));
    curP = curGmm.ComponentProportion(curId);
    
    curLegs{curId} = sprintf( ...
        '%.3f × logN(µ = %.3f, σ = %.3f)', curP, curMean, curStd);
    
    curY = curP .* normpdf(curBinEdges(:), curMean, curStd);
    plot(curAx, curBinEdges(:), curY, 'LineWidth', 2);
end

% Histogram
histogram(curAx, ...
    curLog10Areas, ...
    'BinEdges', curBinEdges, ...
    'Normalization', 'pdf', ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);

% Full GMM
plot(curAx, ...
    curBinEdges(:), curGmm.pdf(curBinEdges(:)), ...
    'LineWidth', 2, 'LineStyle', '--', 'Color', 'black');

curLeg = legend(curAx, curLegs, 'Location', 'SouthOutside');
curLeg.Box = 'off';

curFig.Color = 'white';
curFig.Position(3:4) = [315, 300];
curAx.TickDir = 'out';

xlabel(curAx, 'log_{10}(ASI area [µm²])');
ylabel(curAx, 'PDF estimate');

title(curAx, ...
    {info.filename; info.git_repos{1}.hash; curConfig.title}, ...
    'FontWeight', 'normal', 'FontSize', 10);

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
    fprintf('CV between pairs (mean ± std)\n');
    
    for curPairConfig = curPairConfigs
        curCvs = asiT.area(curPairConfig.synIdPairs);
        curCvs = std(curCvs, 0, 2) ./ mean(curCvs, 2);

        fprintf( ...
            '* %s: %f ± %f\n', ...
            curPairConfig.title, ...
            mean(curCvs), std(curCvs));
    end
    
    fprintf('\n');
    
    fprintf('Areas under and between curves\n');
    fprintf('* Learned fraction: %.1f %%\n', 100 * curLearnedFrac);
    fprintf('* Possibly learned fraction: %.1f %%\n', 100 * curOpenFrac);
    fprintf('* Unlearned fraction: %.1f %%\n', 100 * curUnlearnedFrac);
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
curLimX = [0, 1.5];
curLimY = [-1.5, 0.5];
curImSize = [256, 256];
curMethod = 'kde2d';

% NOTE(amotta): The `curMinMap` and `curMaxMap` matrices have the same size
% as the heatmaps of the CV × log10(avg. ASI) space. They contain the ASI
% areas of the smaller and larger synapses, respectively.
curY = linspace(curLimY(1), curLimY(2), curImSize(1));
curX = linspace(curLimX(1), curLimX(2), curImSize(2));
curX(curX >= sqrt(2)) = nan;

curMinMap = (10 .^ curY(:)) .* (1 - curX / sqrt(2));
curMaxMap = (10 .^ curY(:)) .* (1 + curX / sqrt(2));

curTicksX = linspace(curLimX(1), curLimX(2), 4);
curTicksY = linspace(curLimY(1), curLimY(2), 5);

for curConfig = plotConfigs
    curSaSdConfig = ...
        connectEM.Consistency.buildPairConfigs(asiT, curConfig);
    curSaSdConfig = curSaSdConfig(1);
            
    % Connection types
    curSaSdT = table;
    curSaSdT.synIds = curSaSdConfig.synIdPairs;
    curSaSdT.axonClass = asiT.axonClass(curSaSdT.synIds(:, 1));
    curSaSdT.targetClass = asiT.targetClass(curSaSdT.synIds(:, 1));

    curSaSdT.areas = asiT.area(curSaSdT.synIds);
    curSaSdT.log10Avg = log10(mean(curSaSdT.areas, 2));
    curSaSdT.cv = std(curSaSdT.areas, 0, 2) ./ mean(curSaSdT.areas, 2);

    curSaSdT.mapIdx = [ ...
        discretize(curSaSdT.log10Avg, linspace( ...
            curLimY(1), curLimY(2), curImSize(1) + 1)), ...
        discretize(curSaSdT.cv, linspace( ...
            curLimX(1), curLimX(2), curImSize(2) + 1))];
    curSaSdT.mapIdx = sub2ind( ...
        curImSize, curSaSdT.mapIdx(:, 1), curSaSdT.mapIdx(:, 2));
    
    curCtrlConfigs = struct('synIds', {}, 'title', {}, 'tag', {});
    curCtrlConfigs(1).synIds = curSaSdConfig.synIdPairs(:);
    curCtrlConfigs(1).title = 'SASD';
    
    for curCtrlConfig = curCtrlConfigs
        curKvPairs = { ...
            'xLim', curLimX, 'yLim', curLimY, ...
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
        
        %% Evaluation
        fprintf('Evaluation of\n');
        fprintf('%s\n', curConfigTitle);
        fprintf('%s\n', curCtrlTitle);
        fprintf('\n');
        
        curTitle = cell(2, numel(curRegionProps));
        for curRegionId = 1:numel(curRegionProps)
            curMask = curRegionMask == curRegionId;
            curSaSdFrac = sum(curSaSdMap(curMask));
            curCtrlFrac = sum(curCtrlMap(curMask));
            curDiffFrac = sum(curDiffMap(curMask));
            
            curTitle{1, curRegionId} = ...
                sprintf('%.2f %%', 100 * curSaSdFrac);
            curTitle{2, curRegionId} = ...
                sprintf('%.2f %%', 100 * curDiffFrac);
            
            fprintf('  Region %d\n', curRegionId);
            fprintf([ ...
                '    * Fraction of SASD connections: %.2f %%\n', ...
                '    * Fraction of ctrl. connections: %.2f %%\n', ...
                '    * Delta: %.2f %%\n'], ...
                100 * [curSaSdFrac, curCtrlFrac, curDiffFrac]);
            fprintf('\n');
            
            curMask = curSaSdT.regionId == curRegionId;
           [curOutT, ~, curIds] = unique( ...
               curSaSdT(:, {'axonClass', 'targetClass'}), 'rows');
            curOutT.count = accumarray(curIds, curMask);
            
            curOutT.percentOfRegion = ...
                100 * curOutT.count / sum(curOutT.count);
            curOutT.percentOfRegionPvalue = arrayfun( ...
                @(x, p) 1 - binocdf(x, sum(curOutT.count), p), ...
                curOutT.count - 1, accumarray(curIds, 1) / numel(curIds));
            
            curOutT.prePostTypePercent = ...
                100 * accumarray(curIds, curMask, [], @mean);
            curOutT.prePostTypePvalue = arrayfun( ...
                @(x, n) 1 - binocdf(x, n, curSaSdFrac), ...
                curOutT.count - 1, accumarray(curIds, 1));
            
            disp(curOutT)
        end
        
        curTitle = { ...
            'Significance regions'; ...
            strcat(strjoin(curTitle(1, :), ', '), ' of SASD'); ...
            strcat(strjoin(curTitle(2, :), ', '), ' versus control')};
        title(curPValAx, curTitle, 'FontWeight', 'normal', 'FontSize', 10);
        
        %% Evaluate pairs with small and large synapses
        curFig = figure();
       
        subplot(1, 3, 1);
        curAx = curFig.Children(1);
        
        curX = curMinMap(:);
       [curX, curIds] = sort(curX, 'ascend', 'MissingPlacement', 'last');
        
        hold(curAx, 'on');
        plot(gca, log10(curX), cumsum(curCtrlMap(curIds)));
        plot(gca, log10(curX), cumsum(curSaSdMap(curIds)));
        
        ylabel(curAx, 'Fraction of synapse pairs');
        xlabel(curAx, 'log10(Max ASI area of small synapse [µm²])');
        xlim(curAx, curLimY);
        
        subplot(1, 3, 2);
        curAx = curFig.Children(1);
        
        curX = curMaxMap(:);
       [curX, curIds] = sort(curX, 'descend', 'MissingPlacement', 'last');
       
        hold(curAx, 'on');
        plot(curAx, log10(curX), cumsum(curCtrlMap(curIds)));
        plot(curAx, log10(curX), cumsum(curSaSdMap(curIds)));
        
        xlabel(curAx, 'log10(Min ASI area of large synapse [µm²])');
        xlim(curAx, curLimY);
        
        subplot(1, 3, 3);
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
        
        curAx = curAxes(1);
        curPos = curAx.Position;
        curLeg = {'Random pairs', 'Same-axon same-dendrite pairs'};
        curLeg = legend(curAx, curLeg);
        curLeg.Location = 'SouthOutside';
        curAx.Position = curPos;
        
        annotation( ...
            curFig, 'textbox', [0, 0.9, 1, 0.1], ...
            'String', {curConfigTitle; curCtrlTitle});
        
        curFig.Position(3:4) = [1280, 470];
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
        
        xlim(curAx, [0, 0.05]);
        xlabel(curAx, 'p-value threshold');
        ylabel(curAx, 'Fraction of SASD pairs in region');
        
        curFig.Position(3:4) = [330, 260];
        connectEM.Figure.config(curFig, info);
    end
end
