function detectChiasmataSuper(p, chiParam, agglos, startIdx)
    % set version number
    version = chiParam.version;
    outputDir = chiParam.outputDir;
    
    % decide which function to use
    if isfield(chiParam, 'useSphereClustering') && ...
            chiParam.useSphereClustering
        detectFunc = @connectEM.detectChiasmataSphereClustering;
    else
        detectFunc = @connectEM.detectChiasmata;
    end
    
    % add chiasmata detection parameters to `p`
    chiParam = rmfield(chiParam, {'version', 'outputDir'});
    chiParam = cat(2, fieldnames(chiParam), struct2cell(chiParam));
    chiParam = transpose(chiParam);
    
    p = Util.modifyStruct(p, chiParam{:});
    
    for curIdx = startIdx:500:numel(agglos)
        curOutputDir = fullfile( ...
            outputDir, ...
            sprintf('chiasmataX%s_%d', version, floor(curIdx / 100)), ...
            sprintf('visX%s_%d/', version, curIdx));
        mkdir(curOutputDir);
        
        detectFunc( ...
            p, agglos(curIdx).nodes(:, 1:3), ...
            agglos(curIdx).edges, false, curOutputDir);
    end
end
