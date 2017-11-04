function detectChiasmataSuper(p, chiParam, agglos, startIdx)
    % set version number
    numstr = chiParam.version;
    
    % decide which function to use
    if isfield(chiParam, 'useSphereClustering') && ...
            chiParam.useSphereClustering
        detectFunc = @connectEM.detectChiasmataSphereClustering;
    else
        detectFunc = @connectEM.detectChiasmata;
    end
    
    % add chiasmata detection parameters to `p`
    chiParam = cat(2, fieldnames(chiParam), struct2cell(chiParam));
    chiParam = transpose(chiParam);
    
    p = Util.modifyStruct(p, chiParam{:});
    
    for idx = startIdx:500:numel(agglos)
        outputFolder = fullfile( ...
            chiParam.outputDir, ...
            sprintf('chiasmataX%s_%d', numstr, floor(idx / 100)), ...
            sprintf('visX%s_%d/', numstr, idx));
        mkdir(outputFolder);
        
        detectFunc( ...
            p, agglos(idx).nodes(:, 1:3), ...
            agglos(idx).edges, false, outputFolder);
    end
end
