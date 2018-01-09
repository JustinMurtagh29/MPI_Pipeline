function [flights, info] = loadProject( nmlDir, taskIdFile )
%LOADPROJECT Load all files from a project.
% INPUT nmlDir: string
%           Path to folder containing all nml files of the project.
%       taskIdFile: (Optional) string
%           Task IDs to match the files in the nml folder with the
%           submitted tasks. The order can differ due to the wk task ids
%           (i.e. the first nml file is not necessarily the first task from
%           the task generation list).
%           (Default: the order of the flights is determined by how matlab
%           reads the nml files).
% OUTPUT flights: struct
%           Structure containing the flight path queries. If taskIdFile is
%           supplied then each field of flights has length equal to the
%           number of generated tasks which are ordered in the same way as
%           during task generation. Otherwise the order can be random and
%           there are potentially less tasks then generated.
%        info: struct
%           Miscellaneous information about the flights. Currently contains
%           the fields
%           'missingTasks': Tasks that are missing in nmlDir.
%           'emptyTasks': Tasks that do not have any nodes.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

Util.log('Loading nml files from %s.', nmlDir);
flights = connectEM.Flight.loadFromDirs([], nmlDir);

if ~exist('taskIdFile', 'var') || isempty(taskIdFile)
    return
end

Util.log('Restoring task generation order for flights.');

taskIds = flights.filenamesShort;
[taskIdsGen, startPosGen] = connectEM.Chiasma.Util.loadTaskIds(taskIdFile);

if length(taskIds) < size(taskIdsGen, 1)
    warning('There were less tasks in nmlDir than generated.');
end

[wasFound, idx] = ismember(taskIds, taskIdsGen{:,1});
assert(all(wasFound), 'Not all nml files correspond to generated tasks.');

flights = structfun( ...
    @(x)createSortedFlightMember(x, idx, size(taskIdsGen, 1)), ...
    flights, 'uni', 0);
flights.startPosOrig = startPosGen;

info.missingTasks = find(cellfun(@isempty, flights.filenames));
info.emptyTasks = find(cellfun(@isempty, flights.nodes));

end

function fieldNew = createSortedFlightMember(field, idx, maxIdx)
if ~isempty(field)
    fieldNew = cell(maxIdx, 1);
    fieldNew(idx) = field;
else
    fieldNew = [];
end
end