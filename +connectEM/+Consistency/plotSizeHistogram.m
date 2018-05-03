function fig = plotSizeHistogram(info, synT, plotConfigs, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opt = struct;
    opt.binEdges = linspace(0, 1.2, 61);
    opt = Util.modifyStruct(opt, varargin{:});

    fig = figure();
    fig.Color = 'white';
    fig.Position(3:4) = [820, 475];

    ax = axes(fig);
    hold(ax, 'on');
    ax.TickDir = 'out';

    for curPlotConfig = reshape(plotConfigs, 1, [])
        curSynAreas = synT.area(curPlotConfig.synIds);

        histogram( ...
            ax, curSynAreas, ...
            'BinEdges', opt.binEdges, ...
            'Normalization', 'probability', ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2, ...
            'FaceAlpha', 1);
    end

    xlim(ax, opt.binEdges([1, end]));
    xlabel(ax, 'Synapse area (µm²)');
    ylabel(ax, 'Probability');

    legend( ...
        ax, {plotConfigs.title}, ...
        'Location', 'NorthEast', ...
        'Box', 'off');
    title( ...
       {info.filename, info.git_repos{1}.hash}, ...
        'FontWeight', 'normal', 'FontSize', 10);
end
