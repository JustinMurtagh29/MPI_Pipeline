function fig = plotCouplingHistogram(info, synT, plotConfigs, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opt = struct;
    opt.normalization = 'count';
    opt = Util.modifyStruct(opt, varargin{:});
    
    fig = figure();
    fig.Color = 'white';
    fig.Position(3:4) = [160, 235];

    ax = axes(fig);
    hold(ax, 'on');

    ax.YScale = 'log';
    ax.TickDir = 'out';

    for curPlotConfig = reshape(plotConfigs, 1, [])
        curSynT = synT(curPlotConfig.synIds, :);
        
       [~, ~, curCoupling] = unique(curSynT( ...
            :, {'preAggloId', 'postAggloId'}), 'rows');
        curCoupling = accumarray(curCoupling, 1);

        histogram( ...
            ax, curCoupling, ...
            'Normalization', opt.normalization, ...
            'DisplayStyle', 'stairs');
    end
    
    xMax = cat(2, ax.Children.BinEdges);
    xMax = max(xMax(:));
    
    xlim(ax, [0.5, xMax]);
    xticks(ax, 1:(xMax - 0.5));
    xlabel(ax, 'Synapses per connection');

    switch opt.normalization
        case 'count'
            ylabel(ax, 'Occurences');
        case 'probability'
            ylabel(ax, 'Probability');
    end
    
    curYTicks = log10(ax.YLim);
    curYTicks = [floor(curYTicks(1)), ceil(curYTicks(end))];
    yticks(ax, 10 .^ (curYTicks(1):curYTicks(end)));
    yticklabels(ax, arrayfun( ...
        @num2str, yticks(ax), 'UniformOutput', false));
    
    legend( ...
        ax, {plotConfigs.title}, ...
        'Location', 'SouthOutside');
    connectEM.Figure.config(fig, info);
end
