% A chiasma consisting exclusively of valid flight queries (i.e., no
% comment by HiWis and at least one nodes) might be partitioned in an
% invalid manner. If, for examples, the flight queries indicate overlaps
% A → B and B → 0 (dangling), this marks a contradiction.
%
% In this case we want to requery the exit with a dangling flight path (B)
% plus all the queries incoming onto this exit (A).
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;
rng(0);

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

workingDir = fullfile( ...
    rootDir, 'chiasmataSplitting', '20171009T193744-kmb-on-axons-6c');
axonFile = fullfile( ...
    workingDir, 'outputs', '20171026T171916_results.mat');
outputDir = fullfile(workingDir, 'requeries', 'raw-data');

info = Util.runInfo();

%% load data
data = load(axonFile, 'info', 'summary');
taskFiles = data.info.param.dataFiles;
summary = data.summary;

%% find invalid chiasmata
chiTracings = cat(1, summary.tracings);

% only consider unsolved chiasmata
chiSolved = cat(1, summary.solved);
chiTracings = chiTracings(~chiSolved);
clear chiSolved;

% only consider chiasmata with invalid partitions
chiValid = cellfun(@(t) t.partitionIsValid, chiTracings);
chiTracings = chiTracings(~chiValid);
clear chiInvalid;

%% find tasks to requery
requeryTaskIds = cell(size(chiTracings));
keyboard
for curIdx = 1:numel(chiTracings)
    curTracings = chiTracings{curIdx};
    
    curOverlaps = curTracings.overlaps;
    curOverlaps = cell2mat(curOverlaps')';
    
    % find exit with dangling flights
    curDangIds = curOverlaps(~curOverlaps(:, 2), 1);
    
    % find queries connecting to dangling exits
    curIncomingIds = curOverlaps( ...
        ismember(curOverlaps(:, 2), curDangIds), 1);
    
    % build list of tasks to requery
    curRequeryTaskIds = cat(1, curDangIds, curIncomingIds);
    curRequeryTaskIds = unique(curRequeryTaskIds);
    
    % store result
    requeryTaskIds{curIdx} = curTracings.taskIds(curRequeryTaskIds);
end

% collect all
requeryTaskIds = cat(1, requeryTaskIds{:});

%% look up task definitions
taskFiles = reshape(taskFiles, [], 1);
taskData = cellfun(@(f) ...
    load(f, 'queries', 'taskIds', 'taskDef'), ...
    taskFiles, 'UniformOutput', false);
taskData = Util.concatStructs(1, taskData{:});

% find rows for task IDs
[~, taskRows] = ismember(requeryTaskIds, taskData.taskIds);
assert(all(taskRows > 0));

% scramble order (get multiple HiWis per chiasma)
taskRows = taskRows(randperm(numel(taskRows)));

%% resubmit tasks
curDateStr = datestr(now, 30);
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

out = struct;
out.info = info;
out.queries = taskData.queries(taskRows, :);
out.taskDef = taskData.taskDef(taskRows, :);

saveFile = sprintf('%s_data.mat', curDateStr);
Util.saveStruct(fullfile(outputDir, saveFile), out);

taskParam = struct;
taskParam.dataSet = '2012-09-28_ex145_07x2_ROI2017';
taskParam.taskTypeId = '56d6a7c6140000d81030701e';
taskParam.expDomain = 'queriesMHchiasma';
taskParam.expMinVal = 1;
taskParam.instances = 1;
taskParam.team = 'Connectomics department';
taskParam.project = 'queriesMHchiasma';

taskDefFile = sprintf('%s_flightTasks.txt', curDateStr);
connectEM.exportTaskDefinitions( ...
    taskParam, out.taskDef, fullfile(outputDir, taskDefFile));
