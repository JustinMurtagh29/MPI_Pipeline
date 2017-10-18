function detectChiasmataSuper(startidx, p, useSphereClustering, visualization)
    if nargin < 3
        useSphereClustering = false;
    end
    
    agglos = load(p.inputFile);
    if isfield(agglos,'axons')
        agglos = agglos.axons(agglos.indBigAxons);
    else
        agglos = agglos.dendrites(agglos.indBigDends);
    end
    % set version number
    numstr = p.chiasmataVersion;
    
    % decide which function to use
    if useSphereClustering
        detectFunc = @connectEM.detectChiasmataSphereClustering;
    else
        detectFunc = @connectEM.detectChiasmata;
    end
    
    for idx = startidx : 500 : length(agglos)
        outputFolder = fullfile( ...
            '/tmpscratch/kboerg/chiasmata', ...
            sprintf('chiasmataX%s_%d', numstr, floor(idx / 100)), ...
            sprintf('visX%s_%d/', numstr, idx));
        mkdir(outputFolder);
        
        detectFunc( ...
            p, agglos(idx).nodes(:, 1:3), ...
            agglos(idx).edges, visualization, outputFolder);
    end
end
