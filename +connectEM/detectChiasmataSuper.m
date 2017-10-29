function detectChiasmataSuper(startidx, p, useSphereClustering)
    if nargin < 3
        useSphereClustering = false;
    end
    
    agglos = load(p.inputFile);
    agglos.axons = agglos.axons(agglos.indBigAxons);

    % set version number
    numstr = p.chiasmataVersion;
    
    % decide which function to use
    if useSphereClustering
        detectFunc = @connectEM.detectChiasmataSphereClustering;
    else
        detectFunc = @connectEM.detectChiasmata;
    end
    
    for idx = startidx : 500 : length(agglos.axons)
        outputFolder = fullfile( ...
            '/tmpscratch/kboerg/chiasmata', ...
            sprintf('chiasmataX%s_%d', numstr, floor(idx / 100)), ...
            sprintf('visX%s_%d/', numstr, idx));
        mkdir(outputFolder);
        
        detectFunc( ...
            p, agglos.axons(idx).nodes(:, 1:3), ...
            agglos.axons(idx).edges, false, outputFolder);
    end
end
