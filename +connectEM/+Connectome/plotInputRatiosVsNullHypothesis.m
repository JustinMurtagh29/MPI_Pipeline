% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
param = struct;
param.saveFolder = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connName = 'connectome_axons_18_a_ax_spine_syn_clust';

minSynPost = 10;
info = Util.runInfo();

%% loading data
conn = ...
    connectEM.Connectome.load(param, connName);
axonClasses = unique(conn.axonMeta.axonClass);

%% build class connectome
classConnectome = ...
    connectEM.Connectome.buildClassConnectome( ...
        conn, 'targetClasses', [], 'axonClasses', axonClasses);

classConnectome = transpose(classConnectome);
axonClasses = reshape(axonClasses, 1, []);

%% build dendrite class(es)
dendClasses = struct;
dendClasses(1).ids = find( ...
    conn.denMeta.targetClass ~= 'Somata' ...
  & conn.denMeta.targetClass ~= 'AxonInitialSegment' ...
  & conn.denMeta.synCount >= minSynPost);
dendClasses(1).nullIds = dendClasses(1).ids;
dendClasses(1).title = sprintf( ...
    'Dendrites with ≥ %d synapses (n = %d)', ...
    minSynPost, numel(dendClasses(1).ids));
    
%% plot
for curIdx = 1:numel(dendClasses)
    plotAxonClass( ...
        info, classConnectome, ...
        axonClasses, dendClasses(curIdx));
end

%% plotting
function plotAxonClass(info, classConn, axonClasses, dendClass)
    import connectEM.Specificity.calcExpectedFractionDist;
    import connectEM.Specificity.calcFractionChanceProbs;
    classConn = classConn(dendClass.nullIds, :);
    
    ccId  = find(axonClasses == 'Corticocortical');
    tcId  = find(axonClasses == 'Thalamocortical');
    inhId = find(axonClasses == 'Inhibitory');
    
    % Plotting
    fig = figure();
    fig.Color = 'white';
    fig.Position(3:4) = [840, 440];
    
    binEdges = linspace(0, 1, 21);
    histAxes = cell(0);
    pValAxes = cell(0);
    
    %% Plot exc / (exc + inh)
    obsFrac = classConn(:, [ccId, tcId, inhId]);
    obsFrac = sum(obsFrac(:, 1:2), 2) ./ sum(obsFrac, 2);
    obsFrac(isnan(obsFrac)) = 0;
    
   [expFrac, expCount] = calcExpectedFractionDist( ...
        classConn, [ccId, tcId], [ccId, tcId, inhId]);
   [probLow, probHigh] = calcFractionChanceProbs( ...
        classConn, [ccId, tcId], [ccId, tcId, inhId]);
    
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
    
    xlabel(ax, 'Exc / (Exc + Inh)');
    xlim(ax, binEdges([1, end]));
    
    ax.YAxis.Limits(1) = 10 ^ (-0.1);
    ax.YScale = 'log';
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
        calcExpectedFractionDist(classConn, tcId, [tcId, ccId]);
   [probLow, probHigh] = ...
        calcFractionChanceProbs(classConn, tcId, [tcId, ccId]);
    
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
    
    ax.YAxis.Limits(1) = 10 ^ (-0.1);
    ax.YScale = 'log';
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