function out = splitChiasmataMultiPrepare(p, taskGenDir, taskGenId)
    % Written by
    %   Kevin Boergens <kevin.boergens@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % load results from query generation
    dataFile = sprintf('%s_data.mat', taskGenId);
    data = load(fullfile(taskGenDir, dataFile));
    assert(size(data.queries, 1) == size(data.taskDef, 1));
    
    % load webKNOSSOS task IDs
    taskFile = sprintf('%s_flightTaskIDs.txt', taskGenId);
    taskIds = loadTaskIds(fullfile(taskGenDir, taskFile));
    assert(size(data.queries, 1) == numel(taskIds));
    
    % load flight path NMLs
    nmlDir = sprintf('%s_flightPaths/', taskGenId);
    ff = loadFlightPaths(p, fullfile(taskGenDir, nmlDir));
    
    out = struct;
    out.gitInfo = Util.gitInfo();
    
    out.p = p;
    out.taskGenDir = taskGenDir;
    out.taskGenId = taskGenId;
    out.queries = data.queries;
    out.taskDef = data.taskDef;
    out.taskIds = taskIds;
    out.ff = ff;
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