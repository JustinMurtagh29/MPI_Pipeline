function config(fig, info)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    hasInfo = exist('info', 'var') && not(isempty(info));
    
    set(fig, 'Color', 'white');
    
    axes = findobj(fig, 'type', 'axes');
    set(axes, 'Box', 'off', 'TickDir', 'out');
    
    legends = findobj(fig, 'type', 'legend');
    set(legends, 'Box', 'off');
    
    histograms = findobj(fig, 'type', 'histogram');
    set(histograms, 'DisplayStyle', 'stairs', 'LineWidth', 2);
    
    for ax = reshape(axes, 1, [])
        histograms = findobj(ax, 'type', 'histogram');
        if isempty(histograms); continue; end
        
        binEdges = get(histograms, 'BinEdges');
        binEdges = [ ...
            min(cellfun(@min, binEdges)), ...
            max(cellfun(@max, binEdges))];
        ax.XLim = binEdges;
    end
    
    if isscalar(axes)
        title = axes.Title;
        title.FontSize = 10;
        title.FontWeight = 'normal';
        title.String;
    else
        title = findobj(fig, 'flat', 'type', 'TextBox');
        
        if isempty(title) && hasInfo
            title = annotation(fig, 'textbox', [0, 0.9, 1, 0.1]);
        end
        
        assert(isscalar(title));
        title.EdgeColor = 'none';
        title.HorizontalAlignment = 'center';
    end
        
    if hasInfo
        title.String = [{ ...
            info.filename; ...
            info.git_repos{1}.hash}; ...
            title.String];
    end
end