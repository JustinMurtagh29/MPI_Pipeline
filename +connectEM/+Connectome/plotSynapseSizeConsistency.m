% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');

modeConfigs = struct;
modeConfigs(1).name = 'large';
modeConfigs(1).cvRange = [0, 0.4];
modeConfigs(1).avgLogAsiRange = [-0.14, 0.23];

modeConfigs(2).name = 'small';
modeConfigs(2).cvRange = [0, 0.5];
modeConfigs(2).avgLogAsiRange = [-0.97, -0.6];

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, connFile] = ...
    connectEM.Consistency.loadConnectome(param);

[~, synToSynFile] = fileparts(connFile);
synToSynFile = sprintf('%s_synToSynDists.mat', synToSynFile);
synToSynFile = fullfile(fileparts(connFile), synToSynFile);
synToSyn = load(synToSynFile);

% Loading spine head agglomerates
shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

% Loading augmented graph
graph = Graph.load(rootDir);
graph(~graph.borderIdx, :) = [];

borderAreas = fullfile(rootDir, 'globalBorder.mat');
borderAreas = load(borderAreas, 'borderArea2');
borderAreas = borderAreas.borderArea2;

graph.borderArea = borderAreas(graph.borderIdx);
clear borderAreas;

%% Build axon-spine interface areas
asiT = ...
    connectEM.Connectome.buildAxonSpineInterfaces( ...
        param, graph, shAgglos, conn, syn);
asiT(asiT.type ~= 'PrimarySpine', :) = [];

synT = asiT;
synT.isSpine(:) = true;

%% Prepare data
ccAxonIds = conn.axonMeta.id(conn.axonMeta.axonClass == 'Corticocortical');
tcAxonIds = conn.axonMeta.id(conn.axonMeta.axonClass == 'Thalamocortical');
excAxonIds = union(ccAxonIds, tcAxonIds);

plotConfigs = struct;
plotConfigs(1).synIds = find( ...
    synT.isSpine & ismember(synT.preAggloId, excAxonIds));
plotConfigs(1).title = 'excitatory spine synapses';
plotConfigs(1).tag = 'sp';

plotConfigs(2).synIds = find( ...
    synT.isSpine & ismember(synT.preAggloId, tcAxonIds));
plotConfigs(2).title = 'thalamocortical spine synapses';
plotConfigs(2).tag = 'tc sp';

plotConfigs(3).synIds = find( ...
    synT.isSpine & ismember(synT.preAggloId, ccAxonIds));
plotConfigs(3).title = 'corticocortical spine synapses';
plotConfigs(3).tag = 'cc sp';

%% Report synapse sizes
clear cur*;
fprintf('Synapse sizes\n');
for curConfig = plotConfigs
    curSynT = synT(curConfig.synIds, :);
    curLogAsi = log10(curSynT.area);
    curMeanLog10Asi = mean(curLogAsi);
    curStdLog10Asi = std(curLogAsi);
    
    fprintf( ...
        '* log10(%s ASI): %f ± %f (mean ± std)\n', ...
        curConfig.title, curMeanLog10Asi, curStdLog10Asi)
end

%% Plot distribution of synapse size
connectEM.Consistency.plotSizeHistogram( ...
    info, synT, plotConfigs(1), 'scale', 'log');
connectEM.Consistency.plotSizeHistogram( ...
    info, synT, plotConfigs(2:3), 'scale', 'log');

%% Report number of occurences for degree of coupling
clear cur*;
curPlotConfig = plotConfigs(1);
curSynT = synT(curPlotConfig.synIds, :);
[~, ~, curDegreeOfCoupling] = unique(curSynT( ...
    :, {'preAggloId', 'postAggloId'}), 'rows');
curDegreeOfCoupling = accumarray(curDegreeOfCoupling, 1);

docT = table;
[docT.degreeOfCoupling, ~, curDocCount] = unique(curDegreeOfCoupling);

docT.occurences = accumarray(curDocCount, 1);
docT(docT.degreeOfCoupling == 1, :) = [];

disp(docT)

%% Plot histogram of degree of coupling
connectEM.Consistency.plotCouplingHistogram( ...
    info, synT, plotConfigs(1), 'normalization', 'count');
connectEM.Consistency.plotCouplingHistogram( ...
    info, synT, plotConfigs(2:3), 'normalization', 'probability');

%% Synapse areas vs. degree of coupling
clear cur*;
curPlotCouplings = 1:4;
curConfigs = arrayfun(@(curConfig) struct( ...
    'synT', synT, 'plotConfig', curConfig), plotConfigs);

for curConfig = reshape(curConfigs, 1, [])
    curPlotConfig = curConfig.plotConfig;

    curSynT = curConfig.synT(curPlotConfig.synIds, :);
   [~, ~, curSynT.pairId] = unique(curSynT(:, ...
        {'preAggloId', 'postAggloId'}), 'rows');

    curCouplings = accumarray(curSynT.pairId, 1);
    curSynT.coupling = curCouplings(curSynT.pairId);

    curPlotConfigs = arrayfun( ...
        @(c) struct( ...
            'coupling', c, ...
            'synIds', curPlotConfig.synIds(curSynT.coupling == c), ...
            'title', sprintf('%d-fold %s', c, curPlotConfig.title)), ...
        curPlotCouplings);
    
    connectEM.Consistency.plotSizeHistogram( ...
        info, curConfig.synT, curPlotConfigs, ...
        'scale', 'log', 'title', curPlotConfig.title);

   [curFig, curFit] = ...
        connectEM.Consistency.plotSizeBoxPlot( ...
            info, curConfig.synT, curPlotConfigs, ...
            'title', curPlotConfig.title);
    curFig.Position(3:4) = [250, 420];
    
    curAsiGainVsDoc = 10 ^ curFit.p1;
    fprintf('Fit results for "%s"\n', curPlotConfig.title);
    fprintf('* Fold ASI per degree of coupling: %f\n', curAsiGainVsDoc);
end

%% Illustrate synapse size similarity
clear cur*;
curPlotConfig = plotConfigs(1);
curPairConfigs = ...
    connectEM.Consistency.buildPairConfigs(synT, curPlotConfig);

for curPairConfig = curPairConfigs(1:(end - 1))
    curPairConfig(2) = curPairConfigs(end); %#ok
    curFig = connectEM.Consistency.plotVariabilityPaired( ...
        info, synT, curPlotConfig, curPairConfig, 'lineCount', 10);
    curFig.Position(3:4) = [700, 330];
end

%% Synapse area variability
clear cur*;
curConfigs = struct('synT', synT, 'plotConfigs', plotConfigs);

for curConfig = curConfigs
    for curPlotConfig = curConfig.plotConfigs
        curPairConfigs = ...
            connectEM.Consistency.buildPairConfigs( ...
                curConfig.synT, curPlotConfig);
        curFig = ...
            connectEM.Consistency.plotVariabilityHistogram( ...
                info, curConfig.synT, curPlotConfig, curPairConfigs(:));
        curFig.Position(3:4) = [370, 540];

       [curLearnedFrac, curUnlearnedFrac, curCvThresh] = ...
            connectEM.Consistency.calculateLearnedFraction( ...
                curConfig.synT, curPairConfigs(1), curPairConfigs(end));
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
        connectEM.Consistency.buildPairConfigs(synT, curConfig);
    curSaSdConfig = curPairConfigs(1);
    curRandConfig = curPairConfigs(end);

    curRandSaSdConfig = ...
        connectEM.Consistency.buildPairConfigs( ...
            synT, struct('synIds', curSaSdConfig.synIdPairs(:)));
    curRandSaSdConfig = curRandSaSdConfig(end);
    curRandSaSdConfig.title = sprintf( ...
        'Random pairs from above set (n = %d)', ...
        size(curRandSaSdConfig.synIdPairs, 1));

    for curCtrlConfig = [curRandConfig, curRandSaSdConfig]
        curSaSdT = table;
        curSaSdT.areas = synT.area(curSaSdConfig.synIdPairs);
        curSaSdT.cv = std(curSaSdT.areas, 0, 2) ./ mean(curSaSdT.areas, 2);
        curSaSdT.avgLogAreas = mean(log10(curSaSdT.areas), 2);

        curCtrlT = table;
        curCtrlT.areas = synT.area(curCtrlConfig.synIdPairs);
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
        image(curAx, ...
            uint8(double(intmax('uint8')) ...
          * (1 + curDiffImg / curMaxDiff) / 2));
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

%% Synapse size variability vs. degrees of coupling
clear cur*;
curPlotCouplings = 2:5;
curConfigs = arrayfun(@(curConfig) ...
    struct('synT', synT, 'plotConfig', curConfig), plotConfigs);

for curConfig = curConfigs
    curPlotConfig = curConfig.plotConfig;

    curSynT = curConfig.synT(curPlotConfig.synIds, :);
    curSynT.synId = curPlotConfig.synIds;
    
   [~, ~, curSynT.pairId] = unique(curSynT(:, ...
        {'preAggloId', 'postAggloId'}), 'rows');
    curSynT = sortrows(curSynT, 'pairId');

    curCouplings = accumarray(curSynT.pairId, 1);
    curSynT.coupling = curCouplings(curSynT.pairId);
    
    curPairs = arrayfun( ...
        @(c) struct( ...
            'synIdPairs', transpose(reshape( ...
                curSynT.synId(curSynT.coupling == c), c, [])), ...
            'title', sprintf( ...
                '%d-fold %s', c, curPlotConfig.title)), ...
        curPlotCouplings);
    curCtrlPairs = arrayfun( ...
        @(c) struct( ...
            'synIdPairs', reshape(curPlotConfig.synIds(randperm( ...
                c * floor(numel(curPlotConfig.synIds) / c))), [], c), ...
            'title', sprintf( ...
                '%d-fold %s (random)', c, curPlotConfig.title)), ...
        curPlotCouplings);
    
    curPlotPairs = cat(1, curPairs, curCtrlPairs);
    curPlotConfigs = repelem(curPlotConfig, 1, numel(curPlotCouplings));
    
    curFig = ...
        connectEM.Consistency.plotVariabilityHistogram( ...
            info, curConfig.synT, curPlotConfigs, curPlotPairs);
    curFig.Position(3:4) = [1300, 500];
    
    fprintf('* Significance tests\n');
    fprintf('  p-values for unexpected synapse size similarity\n');
    for curPairConfig = curPlotPairs
        curPValue = ...
            connectEM.Consistency.testVariability( ...
                curConfig.synT, curPairConfig(1), curPairConfig(2));
        fprintf('→ %s: %g\n', curPairConfig(1).title, curPValue);
    end
    
    fprintf('\n');
end
            
%% Variability of largest two synapses
curPlotCouplings = 2:5;
[~, curCouplings, curPlotConfigs] = ...
    ndgrid(1, curPlotCouplings, plotConfigs);

curPairConfigs = @(conf, coup) ...
    connectEM.Consistency.buildLargestPairConfigs(synT, conf, coup);
curPairConfigs = cell2mat(arrayfun( ...
    @(conf, coup) reshape(curPairConfigs(conf, coup), [], 1), ...
    curPlotConfigs, curCouplings, 'UniformOutput', false));

curTitles = arrayfun( ...
    @(coup, conf) sprintf('%d-fold %s', coup, conf.title), ...
    curCouplings, curPlotConfigs, 'UniformOutput', false);
[curPlotConfigs.title] = deal(curTitles{:});

connectEM.Consistency.plotVariabilityHistogram( ...
    info, synT, curPlotConfigs, curPairConfigs);

%% Variability vs. distance
clear cur*;
curConfigs = [ ...
    struct( ...
        'synT', synT, 'plotConfig', plotConfigs(1), ...
        'synToSyn', synToSyn, 'maxDistUm', []), ...
    struct( ...
        'synT', synT, 'plotConfig', plotConfigs(1), ...
        'synToSyn', synToSyn, 'maxDistUm', 20)];
        
for curConfig = curConfigs
    curSynT = curConfig.synT;
    curPlotConfig = curConfig.plotConfig;
    curSynToSyn = curConfig.synToSyn;
    
    curPairConfig = ...
        connectEM.Consistency.buildPairConfigs(curSynT, curPlotConfig);
    curPairConfig = curPairConfig(1);
    
    connectEM.Consistency.plotVariabilityVsDistance( ...
        curSynT, curSynToSyn, curPairConfig, ...
        'maxDistUm', curConfig.maxDistUm);
end
