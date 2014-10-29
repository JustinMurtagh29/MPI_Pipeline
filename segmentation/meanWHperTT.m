%% Read data from BT_tasks
fileID = fopen('/home/mberning/Downloads/BT_tasks.csv');
labels = textscan(fileID, [repmat('%s',1,19) '%*[^\n]'], 1, 'Delimiter', ';');
labels = cellfun(@cell2mat, labels, 'UniformOutput', false);
matrix = textscan(fileID, [repmat('%u',1,12) '%s' repmat('%u',1,3) repmat('%s',1,3) '\n'], 'Delimiter', ';');

%% Read data from BT_tasks
fileID = fopen('/home/mberning/Downloads//BT_uploads.csv');
labels2 = textscan(fileID, [repmat('%s',1,11) '%*[^\n]'], 1, 'Delimiter', ';');
labels2 = cellfun(@cell2mat, labels2, 'UniformOutput', false);
matrix2 = textscan(fileID, [repmat('%u',1,2) '%s%s' repmat('%f',1,3) repmat('%s',1,2) repmat('%u',1,2) '\n'], 'Delimiter', ';');

%% Main loop to determine mean working hours per tasktype
% Look up tasks of tasktype 1 which are finished
tasktypeToCheck = 3;
toKeep = matrix{6} == tasktypeToCheck & matrix{12} == 2;
for i=1:length(matrix)
    matrix3{i} = matrix{i}(toKeep);
end
% Find time spent on last upload
timeSpent = zeros(length(matrix3{1}),1);
for i=1:length(matrix3{1})
    idx = matrix2{2} == matrix3{1}(i);
    if ~(sum(idx) == 0)
        timeSpent(i) = max(matrix2{5}(idx));
    end
end

display(['Tasktype ' num2str(tasktypeToCheck) ': ' num2str(sum(timeSpent == 0)) ' values missing in DB from ' num2str(length(timeSpent)) ' total']);
display(['Tasktype ' num2str(tasktypeToCheck) ': ' num2str(mean(timeSpent(timeSpent ~= 0))) ' hours average for the rest']);