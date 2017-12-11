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
for task = reshape(tasks, 1, [])
    if task.flightDataFileBuild || ~exist(task.flightDataFile, 'file')
        data = connectEM.Chiasma.Heal.loadData( ...
            param, task.genFile, task.idFile, task.nmlDir);
        data.info = info;
        
        Util.saveStruct(task.flightDataFile, data);
        clear data;
    end
end

%%