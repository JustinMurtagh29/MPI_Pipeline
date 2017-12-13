% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;
runId = datestr(now, 30);

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
debugDir = '/home/amotta/Desktop/debug';

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

%% load axons
axonFile = tasks(end).flightDataFile;
axonFile = load(axonFile, 'axonFile');
axonFile = axonFile.axonFile;

axons = load(axonFile, 'axons', 'indBigAxons');
axonIds = find(axons.indBigAxons);
axons = axons.axons(axonIds);

%% determine overlap with axons
flights.overlaps = ...
    connectEM.Flight.overlapWithAgglos( ...
        param, flights, Superagglos.getSegIds(axons), ...
        'minStartEvidence', 13, 'minEndEvidence', 2 * 27);

% override seed agglomerate
[~, queryIdx] = ismember(flights.filenamesShort, queries.taskId);
flights.overlaps(:, 1) = num2cell(queries.aggloId(queryIdx));
flights.seedNodeIds = queries.nodeId(queryIdx);
clear queryIdx;

%% remove invalid flights
% remove flights with comments
[flights, mask] = connectEM.Flight.dropIfCommented(flights);
fprintf('%5.1f %% of flights are commented\n', 100 * mean(~mask));
clear mask;

% remove flights with no or multiple attachments
numOverlaps = cellfun(@numel, flights.overlaps(:, 2));
fprintf('%5.1f %% of flights have no overlap \n', 100 * mean(~numOverlaps));
fprintf('%5.1f %% of flights have too many overlaps \n', 100 * mean(numOverlaps > 1));

%% debugging
if debugDir
    rng(0);
    mkdir(debugDir);
    
    % select random dangling flights
    dangIds = find(~numOverlaps);
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

%%
flights = structfun( ...
    @(f) f(numOverlaps == 1, :), ...
    flights, 'UniformOutput', false);
clear numOverlaps;

flights.overlaps = cell2mat(flights.overlaps(:, 2));

%% assign start agglomerate
% remove duplicate connetions
[~, uniRows] = unique(sort(flights.overlaps, 2), 'rows');
flights = structfun(@(f) f(uniRows, :), flights, 'UniformOutput', false);
clear uniRows;

%% apply
% out = connectEM.Flight.patchIntoAgglos(param, axons, flights);