load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameter.mat', 'p');

temp=load(fullfile(p.saveFolder,'aggloState/axons_06_b.mat'));
temp.axons = temp.axons(temp.indBigAxons);

% get query data, i.e. agglo idx, center idx
queryDir = fullfile(p.saveFolder, 'chiasmata');
queries = connectEM.generateQueriesFromChiasmata(queryDir, temp);

% load task IDs corresponding to queries
taskIdFile = '/tmpscratch/kboerg/IDs_MBKMB_L4_chiasma_axon_queries_26_09_2017.txt';
taskIds = loadTaskIds(taskIdFile);

% sanity check
assert(size(queries, 1) == numel(taskIds));

% get FF structure
nmlDirs = {'/tmpscratch/kboerg/L4_chiasma_axon_queries_26_09_2017_nmls/'};
[ff.segIds, ff.neighbours, ff.filenames, ff.nodes, ff.startNode, ff.comments] = connectEM.lookupNmlMulti(p, nmlDirs, false);
second = @(x)x(2);
ff.filenamesShort = cellfun(@(x)second(strsplit(x,'__')), ff.filenames);

outputDir = '/tmpscratch/kboerg/chiasmaSplit26';

% NOTE(amotta): Sort queries to make processing easier
[~, sortedRows] = sortrows(queries(:, 1:3));
queries = queries(sortedRows, :);
taskIds = taskIds(sortedRows);
clear sortedRows;

for curAxonIdx = 1 : length(temp.axons)
    curAxonIdx
    
    % find all queries in current axon
    curQueryIds = find(queries(:, 1) == curAxonIdx);
    if isempty(curQueryIds); continue; end
    
    curQueries = queries(curQueryIds, :);
    assert(issorted(curQueries(:, 1:3), 'rows'));
    
    % find task IDs and ff indices
    curTaskIds = taskIds(curQueryIds);
    curFfIds = cellfun(@(x) find(strcmp(ff.filenamesShort, x)), curTaskIds);
    assert(numel(curFfIds) == size(curQueries, 1));
    
    % group flights by chiasma
    curFfGroups = accumarray( ...
        curQueries(:, 2), curFfIds, ...
       [curQueries(end, 2), 1], @(ids) {ids}, {});
    curCenterIds = unique(curQueries(:, 7), 'stable');
    
    % create task defintions
    tasks = arrayfun(@(cIdx, ffIds) struct( ...
        'tracings', arrayfun(@(ffIdx) struct( ...
            'segids', ff.segIds(ffIdx), ...
            'nodes', ff.nodes(ffIdx)), ...
            ffIds{1}), ...
        'centeridx', cIdx), ...
        curCenterIds, curFfGroups);
    
    % run splitting
    temp.axons(curAxonIdx).nodesScaled = bsxfun( ...
        @times, temp.axons(curAxonIdx).nodes(:, 1:3), p.raw.voxelSize);
    connectEM.splitChiasmataMulti( ...
        p, temp.axons(curAxonIdx), tasks, ...
        fullfile(outputDir, num2str(curAxonIdx)));
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