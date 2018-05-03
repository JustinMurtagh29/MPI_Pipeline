function fig = plotVariabilityHistogram( ...
        info, synT, plotConfig, pairConfigs, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opt = struct;
    opt.binEdges = linspace(0, 1.5, 21);
    opt = Util.modifyStruct(opt, varargin{:});
    
    fig = figure();
    fig.Color = 'white';
    fig.Position(3:4) = [750, 420];
    
    ax = axes(fig);
    ax.TickDir = 'out';
    axis(ax, 'square');
    hold(ax, 'on');
    
    for curConfig = reshape(pairConfigs, 1, [])
        curCvs = synT.area(curConfig.synIdPairs);
        curCvs = std(curCvs, 0, 2) ./ mean(curCvs, 2);
        
        histogram( ...
            ax, curCvs, ...
            'BinEdges', opt.binEdges, ...
            'Normalization', 'probability', ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);
    end
    
    xlim(ax, opt.binEdges([1, end]));
    xlabel(ax, 'Coefficient of variation');
    ylabel(ax, 'Probability');
    
    legend( ...
        ax, {pairConfigs.title}, ...
        'Location', 'EastOutside', 'Box', 'off');
    title( ...
        ax, {info.filename; info.git_repos{1}.hash; plotConfig.title}, ...
        'FontWeight', 'normal', 'FontSize', 10);
end
