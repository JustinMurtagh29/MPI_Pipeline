function fig = plotSizeBoxPlot(info, synT, plotConfigs, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opt = struct;
    opt.scale = 'log';
    opt.boxWidth = 0.75;
    opt.barWidth = 0.9;
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
    
    %% Preparing data
    coupling = vertcat(plotConfigs.coupling);
    synIds = reshape({plotConfigs.synIds}, [], 1);
    
    dataT = table;
    dataT.coupling = repelem(coupling, cellfun(@numel, synIds));
    dataT.synArea = synT.area(cell2mat(synIds));
    
    % Add jitter for plotting
    dataT.plotX = dataT.coupling ...
        + opt.boxWidth * (rand(height(dataT), 1) - 0.5);
    
    synAreaFit = fit(dataT.coupling, dataT.synArea, 'poly1');

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
    for curConfig = reshape(plotConfigs, 1, [])
        curCoupling = curConfig.coupling;
        
        curX = curCoupling + opt.barWidth * [-0.5, +0.5];
        curY = mean(dataT.synArea(dataT.coupling == curCoupling));
        
        plot( ...
            ax, curX, [curY, curY], ...
            'Color', ax.ColorOrder(2, :), ...
            'LineWidth', 2);
    end
    
    % Interpolation
    plot( ...
        ax, ax.XLim, synAreaFit(ax.XLim), ...
        'Color', 'black', 'LineWidth', 2);
    
    ax.TickDir = 'out';
    ylabel(ax, yLabelText);
    xlabel(ax, 'Degree of coupling');
    
    title( ...
       {info.filename, info.git_repos{1}.hash}, ...
        'FontWeight', 'normal', 'FontSize', 10);
end
