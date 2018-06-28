function [fig, asiFit] = plotSizeBoxPlot(info, synT, plotConfigs, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opt = struct;
    opt.title = {};
    opt.scale = 'log';
    opt.boxWidth = 0.75;
    opt.barWidth = 0.75;
    opt.markerSize = 12;
    opt.markerType = '.';
    opt = Util.modifyStruct(opt, varargin{:});
    
    yLabelText = 'Synapse area (µm²)';
    
    switch opt.scale
        case 'linear'
            % nothing to do
        case 'log'
            synT.area = log10(synT.area);
            yLabelText = sprintf('log_{10}(%s)', yLabelText);
        otherwise
            error('Invalid scale "%s"', opt.scale);
    end
    
    if ~iscell(opt.title); opt.title = {opt.title}; end
    opt.title = cat(1, {info.filename; info.git_repos{1}.hash}, opt.title);
    
    %% Preparing data
    coupling = vertcat(plotConfigs.coupling);
    synIds = reshape({plotConfigs.synIds}, [], 1);
    
    dataT = table;
    dataT.coupling = repelem(coupling, cellfun(@numel, synIds));
    dataT.synArea = synT.area(cell2mat(synIds));
    
    % Add jitter for plotting
    dataT.plotX = dataT.coupling ...
        + opt.boxWidth * (rand(height(dataT), 1) - 0.5);
    
    asiFit = fit(dataT.coupling, dataT.synArea, 'poly1');

    %% Plotting
    fig = figure();
    fig.Color = 'white';

    ax = axes(fig);
    hold(ax, 'on');
    
    scatter( ...
        ax, dataT.plotX, dataT.synArea, ...
        opt.markerSize, opt.markerType);
    
    ax.XLim = prctile(coupling, [0, 100]) + [-0.5, +0.5];
    
    % Means
    dup = @(x) [x, x];
    for curConfig = reshape(plotConfigs, 1, [])
        curCoupling = curConfig.coupling;
        
        curX = curCoupling + opt.barWidth * [-0.5, +0.5];
        curY = dataT.synArea(dataT.coupling == curCoupling);
        curBox = prctile(curY, [0, 25, 50, 75, 100]);
        
        rectangle(ax, 'Position', [ ...
            curX(1), curBox(2), ...
            curX(end) - curX(1), ...
            curBox(4) - curBox(2)]);
        plot(ax, curX, dup(curBox(3)), 'Color', 'black', 'LineWidth', 2);
        
        curWhiskerX = curCoupling + 0.5 * opt.barWidth * [-0.5, +0.5];
        plot(ax, curWhiskerX, dup(curBox(1)), 'Color', 'black');
        plot(ax, curWhiskerX, dup(curBox(end)), 'Color', 'black');
        
        plot( ...
            ax, dup(curCoupling), [curBox(4), curBox(end)], ...
            'Color', 'black', 'LineStyle', '--');
        plot( ...
            ax, dup(curCoupling), [curBox(1), curBox(2)], ...
            'Color', 'black', 'LineStyle', '--');
    end
    
    % Interpolation
    plot( ...
        ax, ax.XLim, asiFit(ax.XLim), ...
        'Color', 'black', 'LineWidth', 2);
    
    ax.TickDir = 'out';
    ylabel(ax, yLabelText);
    xlabel(ax, 'Degree of coupling');
    
    title(ax, opt.title,'FontWeight', 'normal', 'FontSize', 10);
end
