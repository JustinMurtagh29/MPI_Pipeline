function [map, bw] = densityMap(asiAreaPairs, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opts = struct;
    opts.xLim = [0, 1.5];
    opts.yLim = [-1.5, 0.5];
    opts.mapSize = [301, 301];
    opts.bandWidth = [];
    opts = Util.modifyStruct(opts, varargin{:});
    
   [gridY, gridX] = ndgrid( ...
        linspace(opts.yLim(1), opts.yLim(2), opts.mapSize(1)), ...
        linspace(opts.xLim(1), opts.xLim(2), opts.mapSize(2)));
    grid = horzcat(gridY(:), gridX(:));
    
    log10AvgArea = log10(mean(asiAreaPairs, 2));
    cv = std(asiAreaPairs, 0, 2) ./ mean(asiAreaPairs, 2);
    map = horzcat(log10AvgArea, cv);
    
    mask = ...
        opts.yLim(1) <= map(:, 1) & map(:, 1) <= opts.yLim(2) ...
      & opts.xLim(1) <= map(:, 2) & map(:, 2) <= opts.xLim(2);
    support = horzcat(opts.yLim(:), opts.xLim(:));
    
    bwKvPair = {};
    if ~isempty(opts.bandWidth)
        bwKvPair = {'Bandwidth', opts.bandWidth};
    end
    
   [map, ~, bw] = ksdensity( ...
        map(mask, :), grid, ...
        'BoundaryCorrection', 'reflection', ...
        'Support', support, bwKvPair{:});
    map = reshape(map, opts.mapSize);
    map = mean(mask) * map / sum(map(:));
end
