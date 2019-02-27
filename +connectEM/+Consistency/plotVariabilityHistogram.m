function fig = plotVariabilityHistogram( ...
        info, synT, plotConfigs, pairConfigs, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opt = struct;
    opt.binEdges = linspace(0, 1.5, 16);
    opt = Util.modifyStruct(opt, varargin{:});
    
    % Sanity check
    plotConfigSize = size(plotConfigs);
    pairConfigSize = size(pairConfigs);
    
    assert(isequal( ...
        plotConfigSize(2:end), ...
        pairConfigSize(2:end)));
    
    % Homogenize representation
    plotSize = size(plotConfigs);
    plotSize((end + 1):3) = 1;
    
    plotConfigs = reshape( ...
        plotConfigs, [], plotSize(2), plotSize(3));
    pairConfigs = reshape( ...
        pairConfigs, [], plotSize(2), plotSize(3));
    
    fig = figure();
    fig.Color = 'white';
    
    for curRow = 1:plotSize(3)
        for curCol = 1:plotSize(2)
            curPlotConfig = plotConfigs(1, curCol, curRow);
            curPairConfigs = pairConfigs(:, curCol, curRow);
            
            curAx = curCol + (curRow - 1) * plotSize(2);
            curAx = subplot(plotSize(3), plotSize(2), curAx);
            
            curAx.TickDir = 'out';
            axis(curAx, 'square');
            hold(curAx, 'on');
            
            for curPairConfig = reshape(curPairConfigs, 1, [])
                curCvs = synT.area(curPairConfig.synIdPairs);
                curCvs = std(curCvs, 0, 2) ./ mean(curCvs, 2);

                histogram( ...
                    curAx, curCvs, ...
                    'BinEdges', opt.binEdges, ...
                    'Normalization', 'probability', ...
                    'DisplayStyle', 'stairs', ...
                    'LineWidth', 2, ...
                    'FaceAlpha', 1);
            end
            
            xlim(curAx, opt.binEdges([1, end]));
            xlabel(curAx, 'Coefficient of variation');
            ylabel(curAx, 'Probability');
            
            legend( ...
                curAx, {curPairConfigs.title}, ...
                'Location', 'SouthOutside', ...
                'Box', 'off');
            
            title( ...
                curAx, curPlotConfig.title, ...
                'FontWeight', 'normal', ...
                'FontSize', 10);
        end
    end
    
    annotation(fig, ...
        'textbox', [0, 0.9, 1, 0.1], ...
        'String', {info.filename; info.git_repos{1}.hash}, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center');
end
