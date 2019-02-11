% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

modeConfigs = struct(zeros(1, 0));

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
    curLogeAsi = log10(curSynT.area);
    curMeanLogeAsi = mean(curLogeAsi);
    curStdLogeAsi = std(curLogeAsi);
    
    fprintf( ...
        '* log_e(%s ASI): %f ± %f (mean ± std)\n', ...
        curConfig.title, curMeanLogeAsi, curStdLogeAsi)
end


%% Plot distribution of synapse size
for curScale = ["log", "ln"]
    curFig = connectEM.Consistency.plotSizeHistogram( ...
        info, asiT, plotConfigs(1), 'scale', curScale);
    connectEM.Consistency.plotSizeHistogram( ...
        info, asiT, plotConfigs(2:3), 'scale', curScale);
end

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

for curConfig = plotConfigs
    curPairConfigs = ...
        connectEM.Consistency.buildPairConfigs(asiT, curConfig);
    curSaSdConfig = curPairConfigs(1);
    curRandConfig = curPairConfigs(end);
    
    curConfig.title = sprintf( ...
        '%s (n = %d)', curConfig.title, ...
        size(curSaSdConfig.synIdPairs, 1));

    curRandSaSdConfig = ...
        connectEM.Consistency.buildPairConfigs( ...
            asiT, struct('synIds', curSaSdConfig.synIdPairs(:)));
    curRandSaSdConfig = curRandSaSdConfig(end);
    curRandSaSdConfig.title = sprintf( ...
        'Random pairs from above set (n = %d)', ...
        size(curRandSaSdConfig.synIdPairs, 1));
    
    for curCtrlConfig = [curRandConfig, curRandSaSdConfig]
        curSaSdT = table;
        curSaSdT.areas = asiT.area(curSaSdConfig.synIdPairs);
        curSaSdT.cv = std(curSaSdT.areas, 0, 2) ./ mean(curSaSdT.areas, 2);
        curSaSdT.avgLogAreas = mean(log10(curSaSdT.areas), 2);

        curCtrlT = table;
        curCtrlT.areas = asiT.area(curCtrlConfig.synIdPairs);
        curCtrlT.cv = std(curCtrlT.areas, 0, 2) ./ mean(curCtrlT.areas, 2);
        curCtrlT.avgLogAreas = mean(log10(curCtrlT.areas), 2);

        % Density difference map
        curLimX = [0, 1.5];
        curLimY = [-1.5, 0.5];

        curTicksX = linspace(curLimX(1), curLimX(2), 4);
        curTicksY = linspace(curLimY(1), curLimY(2), 5);

        curImSize = [301, 301];

       [curImGridY, curImGridX] = ndgrid( ...
            linspace(curLimY(1), curLimY(2), curImSize(1)), ...
            linspace(curLimX(1), curLimX(2), curImSize(2)));
        curImGrid = cat(2, curImGridY(:), curImGridX(:));

        curSaSdImg = horzcat( ...
            curSaSdT.avgLogAreas, curSaSdT.cv);
        curMask = ...
            curLimY(1) <= curSaSdImg(:, 1) ...
          & curSaSdImg(:, 1) <= curLimY(2) ...
          & curLimX(1) <= curSaSdImg(:, 2) ...
          & curSaSdImg(:, 2) <= curLimX(2);

        curSaSdImg = ksdensity( ...
            curSaSdImg(curMask, :), curImGrid, ...
            'Support', transpose(cat(1, curLimY, curLimX)), ...
            'BoundaryCorrection', 'reflection');
        curSaSdImg = mean(curMask) * curSaSdImg / sum(curSaSdImg(:));
        curSaSdImg = reshape(curSaSdImg, curImSize);

        curCtrlImg = horzcat( ...
            curCtrlT.avgLogAreas, curCtrlT.cv);
        curMask = ...
            curLimY(1) <= curCtrlImg(:, 1) ...
          & curCtrlImg(:, 1) <= curLimY(2) ...
          & curLimX(1) <= curCtrlImg(:, 2) ...
          & curCtrlImg(:, 2) <= curLimX(2);

        curCtrlImg = ksdensity( ...
            curCtrlImg(curMask, :), curImGrid, ...
            'Support', transpose(cat(1, curLimY, curLimX)), ...
            'BoundaryCorrection', 'reflection');
        curCtrlImg = mean(curMask) * curCtrlImg / sum(curCtrlImg(:));
        curCtrlImg = reshape(curCtrlImg, curImSize);

        curMax = max(max(curSaSdImg(:)), max(curCtrlImg(:)));
        curDiffImg = curSaSdImg - curCtrlImg;
        curMaxDiff = max(abs(curDiffImg(:)));

        %% Plotting
        curFig = figure();
        curFig.Color = 'white';
        curFig.Position(3:4) = [360, 1020];

        curAx = subplot(3, 1, 1);
        image(curAx, ...
            uint8(double(intmax('uint8')) ...
          * curSaSdImg / curMax));
        colormap(curAx, jet(256));

        curBar = colorbar('peer', curAx);
        curBar.Ticks = curBar.Limits;
        curBar.TickLabels = {'0', sprintf('%.3g', curMax)};
        curBar.Label.String = {'Fraction of'; curSaSdConfig.title};

        curAx = subplot(3, 1, 2);
        image(curAx, ...
            uint8(double(intmax('uint8')) ...
          * curCtrlImg / curMax));
        colormap(curAx, jet(256));

        curBar = colorbar('peer', curAx);
        curBar.Ticks = curBar.Limits;
        curBar.TickLabels = {'0', sprintf('%.3g', curMax)};
        curBar.Label.String = {'Fraction of'; curCtrlConfig.title};

        curAx = subplot(3, 1, 3);
        hold(curAx, 'on');
        
        curLevels = linspace(-curMaxDiff, +curMaxDiff, 7);
        curLevels((1 + end) / 2) = [];
        
        image(curAx, ...
            uint8(double(intmax('uint8')) ...
          * (1 + curDiffImg / curMaxDiff) / 2));
        contour(curAx, curDiffImg, curLevels, 'LineColor', 'black');
        
        colormap(curAx, jet(256));
        curBar = colorbar('peer', curAx);
        curBar.Label.String = { ...
            'Fraction of observe pairs'; ...
            'relative to null model'};
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

       [~, curTickIdsX] = ismember(curTicksX, ...
           linspace(curLimX(1), curLimX(2), curImSize(2)));
        curTickLabelsX = arrayfun( ...
            @num2str, curTicksX, 'UniformOutput', false);
       [~, curTickIdsY] = ismember(curTicksY, ...
           linspace(curLimY(1), curLimY(2), curImSize(1)));
        curTickLabelsY = arrayfun( ...
            @num2str, curTicksY, 'UniformOutput', false);

        curColors = {'white', 'white', 'black'};
        curImages = {curSaSdImg, curCtrlImg, curDiffImg};
        curTitles = { ...
            curSaSdConfig.title, curCtrlConfig.title, ...
            'Difference in probability density'};

        curAxes = reshape(flip(findobj(curFig, 'Type', 'Axes')), 1, []);
        arrayfun(@(ax) hold(ax, 'on'), curAxes);
        
        for curModeConfig = modeConfigs
            fprintf('# Evaluation of mode "%s"\n', curModeConfig.name);

            curX = curModeConfig.cvRange - curLimX(1);
            curX = curX / (curLimX(end) - curLimX(1));
            curX = round(0.5 + (curImSize(2) - 1) * curX);
            curY = curModeConfig.avgLogAsiRange - curLimY(1);
            curY = curY / (curLimY(end) - curLimY(1));
            curY = round(0.5 + (curImSize(1) - 1) * curY);

            for curIdx = 1:numel(curAxes)
                curAx = curAxes(curIdx);
                curTitle = curTitles{curIdx};
                curColor = curColors{curIdx};

                plot(curAx, ...
                    curX([1, 2, 2, 1, 1]), ...
                    curY([1, 1, 2, 2, 1]), ...
                    'Color', curColor, ...
                    'LineStyle', '--', ...
                    'LineWidth', 2);

                curProb = sum(sum(curImages{curIdx}( ...
                    curY(1):curY(end), curX(1):curX(end))));
                fprintf('* %s: %f %%\n', curTitle, 100 * curProb);
            end

            fprintf('\n');
        end

        set(curAxes, ...
            'Box', 'off', ...
            'TickDir', 'out', ...
            'YDir', 'normal', ...
            'YTick', curTickIdsY, ...
            'YTickLabels', curTickLabelsY, ...
            'XTick', [], ...
            'PlotBoxAspectRatio', [1, 1, 1], ...
            'DataAspectRatioMode', 'auto');
        set(curAxes(end), ...
            'XTick', curTickIdsX, ...
            'XTickLabels', curTickLabelsX);

        arrayfun(@(ax) ylabel(ax, ...
            'Average log_{10}(ASI area [µm²])'), curAxes);
        xlabel(curAxes(end), 'Coefficient of variation');

        annotation( ...
            curFig, ...
            'textbox', [0, 0.9, 1, 0.1], ...
            'String', { ...
                info.filename; info.git_repos{1}.hash; ...
                curConfig.title; curCtrlConfig.title}, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    end
end
