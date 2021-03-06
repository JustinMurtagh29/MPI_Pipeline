function config(fig, info)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    hasInfo = exist('info', 'var') && not(isempty(info));
    
    set(fig, 'Color', 'white');
    
    axes = findobj(fig, 'type', 'axes');
    set(axes, 'Box', 'off', 'TickDir', 'out');
    
    axTitles = cat(1, axes.Title);
    set(axTitles, 'FontWeight', 'normal');
    
    cbars = findobj(fig, 'type', 'colorbar');
    set(cbars, 'Box', 'off', 'TickDir', 'out');
    
    % NOTE(amotta): findobj(fig, 'type', 'numericruler') doesn't work!
    rulers = cat(1, cbars, axes.XAxis, axes.YAxis);
    set(rulers, 'Color', 'black', 'LineWidth', 1);
    
    legends = findobj(fig, 'type', 'legend');
    set(legends, 'Box', 'off');
    
    set( ...
        findobj(fig, 'type', 'histogram'), ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2, ...
        'FaceAlpha', 1);
    
    for ax = reshape(axes, 1, [])
        histograms = findobj(ax, 'type', 'histogram');
        if isempty(histograms); continue; end
        
        binEdges = get(histograms, {'BinEdges'});
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
        annotationPane = findall(fig, 'type', 'AnnotationPane');
        title = findobj(annotationPane, 'type', 'TextBox');
        
        if isempty(title) && hasInfo
            title = annotation(fig, 'textbox', [0, 0.9, 1, 0.1]);
        end
        
        set(title, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    end
        
    if hasInfo
        title.String = [{ ...
            info.filename; ...
            info.git_repos{1}.hash}; ...
            title.String];
    end
    
    fonts = findall(fig, '-property', 'FontName');
    set(cat(1, fonts, rulers), 'FontName', 'Arial', 'FontSize', 10);
end
