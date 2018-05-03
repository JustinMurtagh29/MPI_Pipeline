function fig = plotCouplingHistogram(info, synT, axonClasses, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opt = struct;
    opt.normalization = 'count';
    opt = Util.modifyStruct(opt, varargin{:});
    
    fig = figure();
    fig.Color = 'white';
    fig.Position(3:4) = [495, 400];

    ax = axes(fig);
    axis(ax, 'square');
    hold(ax, 'on');

    ax.YScale = 'log';
    ax.TickDir = 'out';

    for curClassIdx = 1:numel(axonClasses)
        curSynT = axonClasses(curClassIdx).synIds;
        curSynT = synT(curSynT, :);
        
       [~, ~, curCoupling] = unique(curSynT( ...
            :, {'preAggloId', 'postAggloId'}), 'rows');
        curCoupling = accumarray(curCoupling, 1);

        histogram( ...
            ax, curCoupling, ...
            'Normalization', opt.normalization, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2, ...
            'FaceAlpha', 1);
    end

    xMax = max(arrayfun( ...
        @(h) h.BinEdges(end), ax.Children));
    xlim(ax, [0.5, xMax]);
    xticks(ax, 1:(xMax - 0.5));
    xlabel(ax, 'Synapses per connection');

    switch opt.normalization
        case 'count'
            ylabel(ax, 'Occurences');
        case 'probability'
            ylabel(ax, 'Probability');
    end
    
    yticklabels(ax, arrayfun( ...
        @num2str, yticks(ax), ...
        'UniformOutput', false));
    
    legend( ...
        ax, {axonClasses.title}, ...
        'Location', 'NorthEast', 'Box', 'off');
    title( ...
        ax, {info.filename; info.git_repos{1}.hash}, ...
        'FontWeight', 'normal', 'FontSize', 10);
end
