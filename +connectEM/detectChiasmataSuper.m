function detectChiasmataSuper(startidx, p, useSphereClustering, type)
    if nargin < 3
        useSphereClustering = false;
    end
    % Note: For Dendrites set type to true.
    if nargin < 4
        type = false;
    end
    
    if type
        agglos = load(p.inputFile, 'axons', 'indBigAxons');
        agglos = agglos.axons(agglos.indBigAxons);
    else
        agglos = load(p.inputFile, 'dendrites', 'indBigDends');
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
            p.outputDir, ...
            sprintf('chiasmataX%s_%d', numstr, floor(idx / 100)), ...
            sprintf('visX%s_%d/', numstr, idx));
        mkdir(outputFolder);
        
        detectFunc( ...
            p, agglos(idx).nodes(:, 1:3), ...
            agglos(idx).edges, false, outputFolder);
    end
end
