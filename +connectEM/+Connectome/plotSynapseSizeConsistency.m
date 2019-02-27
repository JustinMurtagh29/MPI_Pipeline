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
for curConfig = transpose(plotConfigs(:))
    curSynT = asiT(curConfig.synIds, :);
    
    curLog10MeanAsi = log10(mean(curSynT.area));
    curMeanLog10Asi = mean(log10(curSynT.area));
    curStdLog10Asi = std(log10(curSynT.area));
    
    fprintf( ...
        '* log10(Mean %s ASI [µm²]): %f\n', ...
        curConfig.title, curLog10MeanAsi);
    fprintf( ...
        '* log10(%s ASI [µm²]): %f ± %f (mean ± std.)\n', ...
        curConfig.title, curMeanLog10Asi, curStdLog10Asi);
    fprintf('\n');
end

%% Plot distribution of synapse size
connectEM.Consistency.plotSizeHistogram( ...
    info, asiT, plotConfigs(1, :), 'scale', 'log10');

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

% Show region corresponding to smallest X-th quantile.
curSmall = [0.1, 0.5];

curTicksX = linspace(curLimX(1), curLimX(2), 4);
curTicksY = linspace(curLimY(1), curLimY(2), 5);

for curConfig = plotConfigs
    curSaSdConfig = ...
        connectEM.Consistency.buildPairConfigs(asiT, curConfig);
    curSaSdConfig = curSaSdConfig(1);
    
    curCtrlConfigs = struct('synIds', {}, 'title', {}, 'tag', {});
    curCtrlConfigs(1).synIds = curSaSdConfig.synIdPairs(:);
    curCtrlConfigs(1).title = 'SASD';
    
    curSmallArea = asiT.area(curSaSdConfig.synIdPairs);
    curSmallArea = quantile(curSmallArea(:), curSmall(:));
    
    curSmallY = linspace(curLimX(1), curLimX(2), curImSize(2));
    curSmallY(curSmallY >= sqrt(2)) = nan;
    
    curSmallY = log10(curSmallArea .* sqrt(2) ./ (sqrt(2) - curSmallY));
    curSmallY = curImSize(1) * (curSmallY - curLimY(1)) / diff(curLimY);
    
    for curCtrlConfig = curCtrlConfigs
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
        if ~isnan(curSmall)
            curAxes = findobj(curFig, 'type', 'axes');
            curAxes = flip(transpose(curAxes));
            
            for curAx = curAxes
                curPlot = plot( ...
                    curAx, transpose(curSmallY), 'LineWidth', 2, ...
                    'LineStyle', '--', 'Color', 'white');
            end
            
            curTitle = arrayfun( ...
                @(perc) sprintf('%g-th', perc), ...
                100 * curSmall, 'UniformOutput', false);
            curTitle = sprintf( ...
                'Lower %s percentile lines', strjoin(curTitle, ', '));
            title(curAx, ...
                curTitle, 'FontWeight', 'normal', 'FontSize', 10);
        end
        
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
            
            fprintf('  Region %d\n', curRegionId);
            fprintf([ ...
                '    * Fraction of SASD connections: %.2f %%\n', ...
                '    * Fraction of ctrl. connections: %.2f %%\n', ...
                '    * Delta: %.2f %%\n'], ...
                100 * [curSaSdFrac, curCtrlFrac, curDiffFrac]);
            fprintf('\n');
            
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
    end
end
