function fig = plotSizeHistogram(info, synT, plotConfigs, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opt = struct;
    opt.scale = 'linear';
    opt.binEdges = linspace(0, 1.2, 61);
    opt = Util.modifyStruct(opt, varargin{:});
    
    switch opt.scale
        case 'log10'
            synT.area = log10(synT.area);
        case 'linear'
            % nothing to do
        otherwise
            error('Invalid scale "%s"', opt.scale);
    end

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
