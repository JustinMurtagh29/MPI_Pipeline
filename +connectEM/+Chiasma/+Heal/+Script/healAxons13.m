% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;
runId = datestr(now, 30);

%% configuration
debugDir = ''; % set path to generate debug NMLs
generateRequeries = false; % set to true to generate requeries
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

taskGenDir = fullfile( ...
    rootDir, 'chiasmaSplitHealing', ...
    '20171207T112927-healing-axons-13');
assert(logical(exist(taskGenDir, 'dir')));

taskGenIds = {'20171211T115750'};

tasks = cellfun(@(id) struct( ...
    'genFile', fullfile( ...
        taskGenDir, sprintf('%s_taskGeneration.mat', id)), ...
    'idFile', fullfile( ...
        taskGenDir, sprintf('%s_flightTaskIDs.txt', id)), ...
    'nmlDir', fullfile( ...
        taskGenDir, sprintf('%s_flights', id)), ...
    'flightDataFile', fullfile( ...
        taskGenDir, sprintf('%s_flightData.mat', id)), ...
	'flightDataFileBuild', false), ...
    taskGenIds);

info = Util.runInfo();

%% loading parameters
paramFile = fullfile(rootDir, 'allParameter.mat');
param = load(paramFile, 'p');
param = param.p;

%% load axons
axonFile = tasks(end).flightDataFile;
axonFile = load(axonFile, 'axonFile');
axonFile = axonFile.axonFile;

axons = load(axonFile, 'axons', 'indBigAxons');
axonIdsBig = find(axons.indBigAxons);

allAxons = axons.axons;
axons = axons.axons(axonIdsBig);

%% build data file, if necessary
queries = cell(size(tasks));
flights = cell(size(tasks));

for curIdx = 1:numel(tasks)
    curTask = tasks(curIdx);
    
    if curTask.flightDataFileBuild ...
            || ~exist(curTask.flightDataFile, 'file')
        % need to build data
        curData = connectEM.Chiasma.Heal.loadData( ...
            param, curTask.genFile, curTask.idFile, curTask.nmlDir);
        curData.info = info;
        
        Util.saveStruct(curTask.flightDataFile, curData);
    else
        % load cached file
        curData = load(curTask.flightDataFile);
    end
    
    curQueries = curData.endings;
    curQueries.taskId = curData.taskIds;
    
    queries{curIdx} = curQueries;
    flights{curIdx} = curData.flights;
end

queries = cat(1, queries{:});
flights = Util.concatStructs(1, flights{:});

%% try to find flight for each query
fprintf('# queries handed out: %d\n', size(queries, 1));
[~, uniRows] = unique(queries(:, {'aggloId', 'nodeId'}), 'rows', 'stable');

fprintf('# endings queried: %d\n', numel(uniRows));
queries = queries(uniRows, :);
clear uniRows;

%% remove invalid flights
% remove flights with comments
[flights, mask] = connectEM.Flight.dropIfCommented(flights);
fprintf('%5.1f %% of flights are commented\n', 100 * mean(~mask));
clear mask;

%% determine overlap with axons
flights.overlaps = ...
    connectEM.Flight.overlapWithAgglos( ...
        param, flights, Superagglos.getSegIds(axons), ...
        'minStartEvidence', 13, 'minEndEvidence', 2 * 27);

% remove flights with no or multiple attachments
numOverlaps = cellfun(@numel, flights.overlaps(:, 2));
fprintf('%5.1f %% of flights are dangling\n', 100 * mean(~numOverlaps));
fprintf('%5.1f %% of flights have too many overlaps\n', 100 * mean(numOverlaps > 1));

%% debugging
if logical(debugDir)
    rng(0);
    mkdir(debugDir);
    
    % select random dangling flights
    dangIds = find(numOverlaps == 1);
    dangIds = dangIds(randperm(numel(dangIds), 20));
    axonNames = {'Seed axon', 'Reached axon'};
    
    for curIdx = 1:numel(dangIds)
        curId = dangIds(curIdx);
        curSkel = skeleton();
        
        curTaskId = flights.filenamesShort{curId};
        curFlightName = sprintf('Flight (%s)', curTaskId);
        
        curFlightNodes = flights.nodes{curId};
        curSkel = curSkel.addTree(curFlightName, curFlightNodes);

        curAxonIds = cell2mat(flights.overlaps(curId, :)');
        curAxons = axons(curAxonIds);
        
        for curAxonIdx = 1:numel(curAxonIds)
            curAxonId = curAxonIds(curAxonIdx);
            curAxonName = axonNames{min(curAxonIdx, 2)};
            curAxonName = sprintf('%s (%d)', curAxonName, curAxonId);
            
            curAxon = curAxons(curAxonIdx);
            curSkel = curSkel.addTree( ...
                curAxonName, curAxon.nodes(:, 1:3), curAxon.edges);
        end
        
        curSkel = Skeleton.setParams4Pipeline(curSkel, param);
        curSkelName = sprintf('%d_flight-%s.nml', curIdx, curTaskId);
        curSkel.write(fullfile(debugDir, curSkelName));
    end
end

%% filter flights
flights = structfun( ...
    @(f) f(numOverlaps < 2, :), ...
    flights, 'UniformOutput', false);
clear numOverlaps;

[~, queries.flightId] = ismember(queries.taskId, flights.filenamesShort);
fprintf('# endings answered: %d\n', sum(logical(queries.flightId)));

% fill in seed agglo and node
[~, queryIdx] = ismember(flights.filenamesShort, queries.taskId);
flights.overlaps(:, 1) = num2cell(queries.aggloId(queryIdx));
flights.seedNodeIds = queries.nodeId(queryIdx);
clear queryIdx;

% fill in zero for dangling flights
mask = cellfun(@isempty, flights.overlaps(:, 2));
flights.overlaps(mask, 2) = {0};
clear mask;

flights.overlaps = cell2mat(flights.overlaps);

% self-attachment is forbidden
mask = diff(flights.overlaps, 1, 2) ~= 0;
flights = structfun(@(f) f(mask, :), flights, 'UniformOutput', false);
clear mask;

% drop duplicate connections
[~, uniIds] = unique(sort(flights.overlaps, 2), 'rows');
flights = structfun(@(f) f(uniIds, :), flights, 'UniformOutput', false);
clear uniIds;

%% generate requeries
if generateRequeries
    endings = queries(~queries.flightId, {'aggloId', 'nodeId'});
    
    taskDefFile = sprintf('%s_flightTasks.txt', runId);
    taskDefFile = fullfile(taskGenDir, taskDefFile);

   [taskDefs, endings] = ...
        connectEM.Chiasma.Heal.generateTasks( ...
            param, endings, axons, taskDefFile);
        
	out = struct;
    out.info = info;
    out.endings = endings;
    out.taskDefs = taskDefs;
    out.taskDefFile = taskDefFile;

    outFile = sprintf('%s_taskGeneration.mat', runId);
    Util.saveStruct(fullfile(taskGenDir, outFile), out);
end

%% patching flights into agglomerates
tic;
fprintf('\n');
fprintf('Patching flights into agglomerates... ');
out = connectEM.Flight.patchIntoAgglos(param, allAxons, flights);
fprintf('done!\n');
toc;
