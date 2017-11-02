function queries = loadQueries(nmlDir)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    nmlFiles = NML.findFiles(nmlDir);
    nmlFiles = reshape(nmlFiles, [], 1);

    queries = table;
    queries.nmlFile = fullfile(nmlDir, nmlFiles);

    % parse NML files
   [queries.exits, queriesValid] = cellfun( ...
        @connectEM.Chiasma.Ortho.parseQuery, ...
        queries.nmlFile, 'UniformOutput', false);

    queriesValid = cat(1, queriesValid{:});
    queriesValid = struct2table(queriesValid);
    queries = cat(2, queries, queriesValid);
    
    % extract task IDs
    taskIds = cellfun( ...
        @(s) strsplit(s, '__'), ...
        nmlFiles, 'UniformOutput', false);
    taskIds = cat(1, taskIds{:});
    taskIds = taskIds(:, 2);
    
    queries.taskId = taskIds;
end