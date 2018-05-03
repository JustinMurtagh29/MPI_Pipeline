function fig = plotSynapseSizeHistogram(info, synT, axonClasses, varargin)
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

    for curClassIdx = 1:numel(axonClasses)
        curSynAreas = axonClasses(curClassIdx).synIds;
        curSynAreas = synT.area(curSynAreas);

        histogram( ...
            ax, curSynAreas, ...
            'BinEdges', opt.binEdges, ...
            'Normalization', 'probability', ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2, ...
            'FaceAlpha', 1);
    end

    xlim(ax, binEdges([1, end]));
    xlabel(ax, 'Synapse area (µm²)');
    ylabel(ax, 'Probability');

    legend( ...
        ax, {axonClasses.title}, ...
        'Location', 'NorthEast', ...
        'Box', 'off');
    title( ...
       {info.filename, info.git_repos{1}.hash}, ...
        'FontWeight', 'normal', 'FontSize', 10);
end
