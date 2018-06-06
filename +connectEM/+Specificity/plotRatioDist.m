function plotRatioDist(info, classConn, axonClasses, dendClass)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    import connectEM.Specificity.calcExpectedRatioDist;
    import connectEM.Specificity.calcRatioChanceProbs;
    classConn = classConn(dendClass.nullIds, :);
    
    ccId  = find(axonClasses == 'Corticocortical');
    tcId  = find(axonClasses == 'Thalamocortical');
    inhId = find(axonClasses == 'Inhibitory');
    
    % Plotting
    fig = figure();
    fig.Color = 'white';
    fig.Position(3:4) = [840, 440];
    
    binEdges = linspace(0, 1, 51);
    histAxes = cell(0);
    pValAxes = cell(0);
    
    %% Plot inh / (exc + inh)
    obsFrac = classConn(:, [ccId, tcId, inhId]);
    obsFrac = sum(obsFrac(:, 3), 2) ./ sum(obsFrac, 2);
    obsFrac(isnan(obsFrac)) = 0;
    
   [expFrac, expCount] = calcExpectedRatioDist( ...
        classConn, inhId, [ccId, tcId, inhId]);
   [probLow, probHigh] = calcRatioChanceProbs( ...
        classConn, inhId, [ccId, tcId, inhId]);
    
    binId = discretize(expFrac, binEdges);
    binCount = accumarray(binId, expCount);
    
    ax = subplot(2, 2, 1);
    axis(ax, 'square');
    hold(ax, 'on');
    
    histogram(ax, ...
        obsFrac, ...
        'BinEdges', binEdges, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
    histogram(ax, ...
        'BinEdges', binEdges, ...
        'BinCounts', binCount, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
    
    xlabel(ax, 'Inh / (Exc + Inh)');
    xlim(ax, binEdges([1, end]));
    
    ax.TickDir = 'out';
    histAxes{end + 1} = ax;
    
    ax = subplot(2, 2, 3);
    colorOrder = ax.ColorOrder;
    axis(ax, 'square');
    
    yyaxis(ax, 'right');
    hold(ax, 'on');
    
    probLow = sort(probLow, 'ascend');
    probLowRatio = (1:numel(probLow)) ./ numel(probLow);
    probLowRatio = probLowRatio(:) ./ probLow;
    
    probHigh = sort(probHigh, 'ascend');
    probHighRatio = (1:numel(probHigh)) ./ numel(probHigh);
    probHighRatio = probHighRatio(:) ./ probHigh;
    
    plot(ax, ...
        probLow, probLowRatio, ...
        'Color', colorOrder(1, :), ...
        'LineStyle', '-');
    plot(ax, ...
        probHigh, probHighRatio, ...
        'Color', colorOrder(2, :), ...
        'LineStyle', '-');
    plot(ax, [0, 1], [1, 1], 'k--');
    
    xlim(ax, [0, 1]);
    ylim(ax, [0, 2]);
    
    yyaxis(ax, 'left');
    hold(ax, 'on');
    
    histogram(ax, ...
        probLow, ...
        'BinEdges', binEdges, ...
        'DisplayStyle', 'stairs', ...
        'EdgeColor', colorOrder(1, :), ...
        'LineWidth', 2);
    histogram(ax, ...
        probHigh, ...
        'BinEdges', binEdges, ...
        'DisplayStyle', 'stairs', ...
        'EdgeColor', colorOrder(2, :), ...
        'LineWidth', 2);
    xlabel(ax, 'Probability');
    
    ax.TickDir = 'out';
    pValAxes{end + 1} = ax;
    
    %% Plot expected tc / (tc + cc)
    obsFrac = classConn(:, [tcId, ccId]);
    obsFrac = obsFrac(:, 1) ./ sum(obsFrac, 2);
    obsFrac(isnan(obsFrac)) = 0;
    
   [expFrac, expCount] = ...
        calcExpectedRatioDist(classConn, tcId, [tcId, ccId]);
   [probLow, probHigh] = ...
        calcRatioChanceProbs(classConn, tcId, [tcId, ccId]);
    
    binId = discretize(expFrac, binEdges);
    binCount = accumarray(binId, expCount);
    
    ax = subplot(2, 2, 2);
    axis(ax, 'square');
    hold(ax, 'on');
    
    histogram(ax, ...
        obsFrac, ...
        'BinEdges', binEdges, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
    histogram(ax, ...
        'BinEdges', binEdges, ...
        'BinCounts', binCount, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
    
    xlabel(ax, 'TC / (TC + CC)');
    xlim(ax, binEdges([1, end]));
    
    ax.TickDir = 'out';
    histAxes{end + 1} = ax;
    
    % Legend
    axPos = ax.Position;
    legend(ax, ...
        'Observed', ...
        'Expected (multinomial)', ...
        'Location', 'North');
    ax.Position = axPos;
    
    ax = subplot(2, 2, 4);
    colorOrder = ax.ColorOrder;
    axis(ax, 'square');
    
    yyaxis(ax, 'right');
    hold(ax, 'on');
    
    probLow = sort(probLow, 'ascend');
    probLowRatio = (1:numel(probLow)) ./ numel(probLow);
    probLowRatio = probLowRatio(:) ./ probLow;
    
    probHigh = sort(probHigh, 'ascend');
    probHighRatio = (1:numel(probHigh)) ./ numel(probHigh);
    probHighRatio = probHighRatio(:) ./ probHigh;
    
    plot(ax, ...
        probLow, probLowRatio, ...
        'Color', colorOrder(1, :), ...
        'LineStyle', '-');
    plot(ax, ...
        probHigh, probHighRatio, ...
        'Color', colorOrder(2, :), ...
        'LineStyle', '-');
    plot(ax, [0, 1], [1, 1], 'k--');
    
    xlim(ax, [0, 1]);
    ylim(ax, [0, 2]);
    
    yyaxis(ax, 'left');
    hold(ax, 'on');
    
    histogram(ax, ...
        probLow, ...
        'BinEdges', binEdges, ...
        'DisplayStyle', 'stairs', ...
        'EdgeColor', colorOrder(1, :), ...
        'LineWidth', 2);
    histogram(ax, ...
        probHigh, ...
        'BinEdges', binEdges, ...
        'DisplayStyle', 'stairs', ...
        'EdgeColor', colorOrder(2, :), ...
        'LineWidth', 2);
    
    xlabel('Probability');
    ax.TickDir = 'out';
    pValAxes{end + 1} = ax;
    
    % Legend
    axPos = ax.Position;
    legend(ax, ...
        'Prob[F ≤ f]', ...
        'Prob[F ≥ f]', ...
        'Location', 'North');
    
    xlabel(ax, 'Probability');
    ax.Position = axPos;

    annotation(fig, ...
        'textbox', [0, 0.9, 1, 0.1], ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'String', { ...
            info.filename; ...
            info.git_repos{1}.hash; ...
            dendClass.title});
        
    %% Fix axes
    histAxes = cat(1, histAxes{:});
    yLims = sort(cat(2, histAxes.YLim));
   [histAxes.YLim] = deal(yLims([1, end]));
   
    pValAxes = cat(1, pValAxes{:});
    yLims = sort(cat(2, pValAxes.YLim));
   [pValAxes.YLim] = deal(yLims([1, end]));
end
