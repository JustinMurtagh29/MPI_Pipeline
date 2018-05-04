function fig = plotVariabilityPaired( ...
        info, synT, plotConfig, pairConfigs, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opt = struct;
    opt.lineCount = 10;
    opt = Util.modifyStruct(opt, varargin{:});
    
    % Calculate and plot in log-space
    synT.area = log10(synT.area);
    
    % Select samples to plot
   [randAreaPairs, meanPairs] = arrayfun( ...
        @(pairConfig) forPairConfig(opt, synT, pairConfig), ...
        pairConfigs, 'UniformOutput', false);
    
    randAreaPairs = cell2mat(reshape(randAreaPairs, 1, 1, []));
    randAreaPairs = permute(randAreaPairs, [2, 3, 1]);
    randAreaPairs = reshape(randAreaPairs, 2, []);
    randAreaPairs = transpose(randAreaPairs);
    meanPairs = cell2mat(meanPairs(:));
   
    fig = figure();
    fig.Color = 'white';
    
    ax = axes(fig);
    axis(ax, 'square');
    hold(ax, 'on');
    
    plot( ...
        ax, [1, 2], randAreaPairs, ...
        'LineStyle', '--', ...
        'MarkerSize', 9, ...
        'Marker', '.');
    plot( ...
        ax, [1, 1.5, 2], meanPairs, ...
        'LineWidth', 2);
    
    colors = ax.ColorOrder(1:numel(pairConfigs), :);
    colors = num2cell(repmat(colors, 1 + opt.lineCount, 1), 2);
    
    lines = flip(ax.Children);
   [lines.Color] = deal(colors{:});
    
    ax.TickDir = 'out';
    xlim(ax, [0.9, 2.1]);
    xticks(ax, [1, 2]);
    xlabel(ax, 'Synapse');
    ylabel(ax, 'log_{10}(ASI)');
    
    legend( ...
        flip(ax.Children(1:2)), ...
        {pairConfigs.title}, ...
        'Location', 'EastOutside', ...
        'Box', 'off');
    
    title( ...
        ax, {info.filename; info.git_repos{1}.hash; plotConfig.title}, ...
        'FontWeight', 'normal', 'FontSize', 10);
end

function [areaPairs, meanPair] = forPairConfig(opt, synT, pairConfig)
    synIdPairs = pairConfig.synIdPairs;
    
    % Select random examples for plotting
    rng(0);
    areaPairs = randperm(size(synIdPairs, 1), opt.lineCount);
    areaPairs = synT.area(synIdPairs(areaPairs, :));
    
    % Calculate mean pair
    pairMean = mean(synT.area(synIdPairs), 2);
    pairDelta = mean(synT.area(synIdPairs(:, 1)) - pairMean);
    
    % Build output
    meanPair = mean(pairDelta) * [1, 0, -1];
    meanPair = mean(pairMean) + meanPair;
end
