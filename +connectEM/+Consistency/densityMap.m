function [map, bw] = densityMap(asiAreaPairs, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opts = struct;
    opts.xAxis = 'cv';
    opts.xLim = [];
    opts.yScale = 'log10';
    opts.yLim = [];
    opts.mapSize = [301, 301];
    opts.bandWidth = [];
    opts.method = 'ksdensity';
    
    opts = Util.modifyStruct(opts, varargin{:});
    
    opts.xAxis = lower(opts.xAxis);
    opts.yScale = lower(opts.yScale);
    opts.method = lower(opts.method);
    
    switch opts.xAxis
        case 'cv'
            x = std(asiAreaPairs, 0, 2) ./ mean(asiAreaPairs, 2);
            if isempty(opts.xLim); opts.xLim = [0, 1.5]; end
        case 'reldiff'
            x = abs(diff(asiAreaPairs, 1, 2)) ./ mean(asiAreaPairs, 2);
            if isempty(opts.xLim); opts.xLim = [0, 2]; end
        otherwise
            error('Invalid X axis "%s"', opts.xAxis);
    end
    
    switch opts.yScale
        case 'linear'
            y = mean(asiAreaPairs, 2);
            if isempty(opts.yLim); opts.yLim = [0, 1]; end
        case 'log10'
            y = log10(mean(asiAreaPairs, 2));
            if isempty(opts.yLim); opts.yLim = [-1.5, 0.5]; end
        otherwise
            error('Invalid Y scale "%s"', opts.yScale);
    end
    
    map = horzcat(y, x);
    
    mask = ...
        opts.yLim(1) <= map(:, 1) & map(:, 1) <= opts.yLim(2) ...
      & opts.xLim(1) <= map(:, 2) & map(:, 2) <= opts.xLim(2);
    support = horzcat(opts.yLim(:), opts.xLim(:));
    
    switch opts.method
        case 'ksdensity'
           [gridY, gridX] = ndgrid( ...
                linspace(opts.yLim(1), opts.yLim(2), opts.mapSize(1)), ...
                linspace(opts.xLim(1), opts.xLim(2), opts.mapSize(2)));
            grid = horzcat(gridY(:), gridX(:));

            bwKvPair = {};
            if ~isempty(opts.bandWidth)
                bwKvPair = {'Bandwidth', opts.bandWidth};
            end

           [map, ~, bw] = ksdensity( ...
                map(mask, :), grid, ...
                'BoundaryCorrection', 'reflection', ...
                'Support', support, bwKvPair{:});
            map = reshape(map, opts.mapSize);
            
        case 'kde2d'
           [bw, map] = connectEM.Libs.kde2d( ...
                map(mask, :), opts.mapSize(1), ...
                support(1, :), support(2, :), ...
                opts.bandWidth);
            map = transpose(map);
            
        otherwise
            error('Invalid method "%s"', opts.method);
    end   
    
    map = mean(mask) * map / sum(map(:));
end
