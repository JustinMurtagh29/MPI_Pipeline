function [queries, taskIds, ff] = ...
        splitChiasmataMultiPrepare(p, taskGenDir, taskGenId)
    % Written by
    %   Kevin Boergens <kevin.boergens@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % load results from query generation
    queriesFile = sprintf('%s_data.mat', taskGenId);
    queries = load(fullfile(taskGenDir, queriesFile));
    queries = queries.queries;
    
    % load webKNOSSOS task IDs
    taskFile = sprintf('%s_flightTaskIDs.txt', taskGenId);
    taskIds = loadTaskIds(fullfile(taskGenDir, taskFile));
    assert(size(queries, 1) == numel(taskIds));
    
    % load flight path NMLs
    nmlDir = sprintf('%s_flightPaths/', taskGenId);
    ff = loadFlightPaths(p, nmlDir);
end

function taskIds = loadTaskIds(taskFile)
    fid = fopen(taskFile, 'rt');
    data = fread(fid, 'char=>char');
    fclose(fid);
    
    % split into lines
    data = reshape(data, 1, []);
    data = strsplit(data, '\n');
    
    % remove header
    data(1) = [];
    
    % extract first column with task IDs
    taskIds = cellfun( ...
        @(x) x(1:(find(x == ',', 1) - 1)), ...
        data(:), 'UniformOutput', false);
end

function ff = loadFlightPaths(p, nmlDirs)
    if ~iscell(nmlDirs)
        nmlDirs = {nmlDirs};
    end
    
    ff = struct;
   [ff.segIds, ff.neighbours, ff.filenames, ...
    ff.nodes, ff.startNode, ff.comments] = ...
        connectEM.lookupNmlMulti(p, nmlDirs, false);
    
    second = @(x) x(2);
    ff.filenamesShort = cellfun( ...
        @(x) second(strsplit(x, '__')), ff.filenames);
end