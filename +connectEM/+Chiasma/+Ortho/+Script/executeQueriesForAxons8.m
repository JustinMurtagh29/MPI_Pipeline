% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
nmlDir = '/home/amotta/Desktop/annotationZips4468876378594053618L4_chiasma_correction_31_10_2017_nmls';

%% find NML files
nmlFiles = NML.findFiles(nmlDir);
nmlFiles = reshape(nmlFiles, [], 1);

%% parse NML files
queries = table;
queries.nmlFile = fullfile(nmlDir, nmlFiles);

% extract task IDs
queryData = cellfun( ...
    @(s) strsplit(s, '__'), ...
    nmlFiles, 'UniformOutput', false);
queryData = cat(1, queryData{:});
queries.taskId = queryData(:, 2);
clear queryData;

% TODO(amotta):
% extract axon and node IDs based on `taskId`

% parse NML files
[queries.exits, queries.comments] = cellfun( ...
    @connectEM.Chiasma.Ortho.parseQuery, ...
    queries.nmlFile, 'UniformOutput', false);

%% statistics
fprintf('\n');
fprintf('# answered queries: %d\n', size(queries, 1));
fprintf('# commented queries: %d\n', sum( ...
    not(cellfun(@isempty, queries.comments))));