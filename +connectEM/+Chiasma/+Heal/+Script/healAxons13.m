% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;
runId = datestr(now, 30);

%% configuration
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
    
%% remove invalid flights
% remove flights with comments
flights = connectEM.Flight.dropIfCommented(flights);

% remove flights with no or multiple attachments
mask = cellfun(@numel, flights.overlaps(:, 2)) == 1;
flights = structfun(@(f) f(mask, :), flights, 'UniformOutput', false);
clear mask;

flights.overlaps = cell2mat(flights.overlaps(:, 2));