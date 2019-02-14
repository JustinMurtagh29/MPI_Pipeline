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

plotConfigs = struct;
plotConfigs(1).synIds = find( ...
    asiT.type == 'PrimarySpine' & ismember( ...
    asiT.axonClass, {'Corticocortical', 'Thalamocortical'}));
plotConfigs(1).title = 'excitatory primary spine synapses';
plotConfigs(1).tag = 'exc pri sp';

plotConfigs(2).synIds = find( ...
    asiT.type == 'PrimarySpine' ...
  & asiT.axonClass == 'Thalamocortical');
plotConfigs(2).title = 'thalamocortical primary spine synapses';
plotConfigs(2).tag = 'tc pri sp';

plotConfigs(3).synIds = find( ...
    asiT.type == 'PrimarySpine' ...
  & asiT.axonClass == 'Corticocortical');
plotConfigs(3).title = 'corticocortical primary spine synapses';
plotConfigs(3).tag = 'cc pri sp';

%% Report synapse sizes
clear cur*;
fprintf('\nSynapse sizes\n');
for curConfig = plotConfigs
    curSynT = asiT(curConfig.synIds, :);
    curLog10MeanAsi = log10(mean(curSynT.area));
    
    fprintf( ...
        '* log_10(Mean %s ASI [µm²]): %f\n', ...
        curConfig.title, curLog10MeanAsi)
end

fprintf('\n');

%% Plot distribution of synapse size
connectEM.Consistency.plotSizeHistogram( ...
    info, asiT, plotConfigs(1), 'scale', 'log10');
connectEM.Consistency.plotSizeHistogram( ...
    info, asiT, plotConfigs(2:3), 'scale', 'log10');

%% Fit mixture of Gaussians to size distribution
clear cur*;
rng(0);

curN = 2;
curConfig = plotConfigs(1);
curBinEdges = linspace(-2, 0.5, 51);

curLog10Areas = log10(asiT.area(curConfig.synIds));
curGmm = fitgmdist(curLog10Areas, curN);

curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');

% GMM components
for curId = 1:curN
    curMean = curGmm.mu(curId);
    curStd = sqrt(curGmm.Sigma(curId));
    curP = curGmm.ComponentProportion(curId);
    
    fprintf('GMM component %d\n', curId);
    fprintf('* Mean: %g\n', curMean);
    fprintf('* Standard deviation: %g\n', curStd);
    fprintf('* Mixing coefficient: %g\n', curP);
    fprintf('\n');
    
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

curFig.Color = 'white';
curFig.Position(3:4) = [315, 225];
curAx.TickDir = 'out';

xlabel(curAx, 'log_{10}(ASI area [µm²])');
ylabel(curAx, 'PDF estimate');

title(curAx, ...
    {info.filename; info.git_repos{1}.hash; curConfig.title}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Report number of occurences for degree of coupling
clear cur*;
curPlotConfig = plotConfigs(1);
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
    info, asiT, plotConfigs(1), 'normalization', 'count');

%% Synapse area variability
clear cur*;
curConfigs = struct('synT', asiT, 'plotConfigs', plotConfigs);

for curConfig = curConfigs
    for curPlotConfig = curConfig.plotConfigs(1)
        curPairConfigs = ...
            connectEM.Consistency.buildPairConfigs( ...
                curConfig.synT, curPlotConfig);
        curFig = ...
            connectEM.Consistency.plotVariabilityHistogram( ...
                info, curConfig.synT, curPlotConfig, curPairConfigs(:));
        curFig.Position(3:4) = [370, 540];

       [curLearnedFrac, curUnlearnedFrac, curCvThresh] = ...
            connectEM.Consistency.calculateLearnedFraction( ...
                curConfig.synT, curPairConfigs(1), curPairConfigs(end), ...
                'method', 'maxdifference');
        curOpenFrac = 1 - curLearnedFrac - curUnlearnedFrac;
        
        % Quantitative reporting
        fprintf('**%s**\n', curPlotConfig.title);
        curSynIdPairs = curPairConfigs(1).synIdPairs;
        
        curSynT = curConfig.synT( ...
            curSynIdPairs(:, 1), {'preAggloId', 'postAggloId'});
        curSynT.areas = curConfig.synT.area(curSynIdPairs);
        curSynT.cv = std(curSynT.areas, 0, 2) ./ mean(curSynT.areas, 2);
        curSynT.low = curSynT.cv < curCvThresh;
        
        curAxonT = curSynT(:, {'preAggloId'});
       [curAxonT, ~, curGroupIds] = unique(curAxonT);
        curAxonT.allConn = accumarray(curGroupIds, 1);
        curAxonT.lowConn = accumarray(curGroupIds, curSynT.low);
        curAxonT = sortrows(curAxonT, 'lowConn', 'descend');
        
        curLowThresh = round(height(curSynT) * curLearnedFrac);
        curAxonFracLow = find(cumsum(curAxonT.lowConn) > curLowThresh, 1);
        curAxonFracLow = curAxonFracLow / height(curAxonT) %#ok
        curAxonFracUp = curLowThresh / height(curAxonT) %#ok
        
        fprintf('CV between pairs (mean ± std)\n');
        for curPairConfig = curPairConfigs
            curCvs = curConfig.synT.area(curPairConfig.synIdPairs);
            curCvs = std(curCvs, 0, 2) ./ mean(curCvs, 2);
            
            fprintf( ...
                '* %s: %f ± %f\n', ...
                curPairConfig.title, ...
                mean(curCvs), std(curCvs));
        end
            
        fprintf('\n');
        
        fprintf('* Learned fraction: %.1f %%\n', 100 * curLearnedFrac);
        fprintf('* Possibly learned fraction: %.1f %%\n', 100 * curOpenFrac);
        fprintf('* Unlearned fraction: %.1f %%\n', 100 * curUnlearnedFrac);
        fprintf('* CV threshold: %.2f\n', curCvThresh);
        fprintf('\n');

        fprintf('Significance tests\n');
        fprintf('p-values for unexpected synapse size similarity\n');
        for curPairConfig = curPairConfigs(1:(end - 1))
            curPValue = ...
                connectEM.Consistency.testVariability( ...
                    curConfig.synT, curPairConfig, curPairConfigs(end));
            fprintf('* %s: %g\n', curPairConfig.title, curPValue);
        end
        
        fprintf('\n');
    end
end

%% Correlation of synapse size correlation with synapse size
clear cur*;

% Density difference map
curLimX = [0, 1.5];
curLimY = [-1.5, 0.5];
curImSize = [256, 256];
curMethod = 'kde2d';

curTicksX = linspace(curLimX(1), curLimX(2), 4);
curTicksY = linspace(curLimY(1), curLimY(2), 5);

for curConfig = plotConfigs
    curSaSdConfig = ...
        connectEM.Consistency.buildPairConfigs(asiT, curConfig);
    curSaSdConfig = curSaSdConfig(1);
    
    curCtrlConfigs = struct('synIds', {}, 'title', {}, 'tag', {});
    curCtrlConfigs(1).synIds = curConfig.synIds(:);
    curCtrlConfigs(1).title = 'all';
    curCtrlConfigs(2).synIds = curSaSdConfig.synIdPairs(:);
    curCtrlConfigs(2).title = 'SASD';
    
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
        end
    end
end
