function fig = plotSizeHistogram(info, synT, plotConfigs, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opt = struct;
    opt.title = {};
    opt.scale = 'linear';
    opt.binEdges = [];
    opt = Util.modifyStruct(opt, varargin{:});
    
    xLabelText = 'Synapse area (µm²)';
    
    switch opt.scale
        case 'linear'
            defaultBinEdges = linspace(0, 1.2, 61);
        case {'ln', 'loge'}
            synT.area = log(synT.area);
            defaultBinEdges = linspace(-4, 1, 51);
            xLabelText = sprintf('log_{e}(%s)', xLabelText);
        case {'log', 'log10'}
            synT.area = log10(synT.area);
            defaultBinEdges = linspace(-2, 0.5, 51);
            xLabelText = sprintf('log_{10}(%s)', xLabelText);
        otherwise
            error('Invalid scale "%s"', opt.scale);
    end
    if isempty(opt.binEdges)
        opt.binEdges = defaultBinEdges;
    end
    
    if ~iscell(opt.title); opt.title = {opt.title}; end
    opt.title = cat(1, {info.filename; info.git_repos{1}.hash}, opt.title);

    fig = figure();
    fig.Color = 'white';
    fig.Position(3:4) = [360, 230];

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
    xlabel(ax, xLabelText);
    ylabel(ax, 'Probability');

    legend(ax, {plotConfigs.title}, 'Location', 'NorthEast', 'Box', 'off');
    title(opt.title, 'FontWeight', 'normal', 'FontSize', 10);
end
