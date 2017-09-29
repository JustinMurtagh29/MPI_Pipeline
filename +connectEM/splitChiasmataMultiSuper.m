function splitChiasmataMultiSuper(p)
    % Written by
    %   Kevin Boergens <kevin.boergens@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    oldAgglos = load(fullfile( ...
        p.saveFolder, 'aggloState', 'axons_06_b.mat'));
    
    agglos = oldAgglos;
    agglos.axons = agglos.axons(agglos.indBigAxons);
    
    curDir = fullfile( ...
        p.saveFolder, 'chiasmataSplitting', ...
        '20170928T132211-projects-a-and-b');
    outputDir = fullfile(curDir, 'outputs');

    load(fullfile(curDir, 'input-data.mat'));
    clear curDir;
    % load redundancy definition
    load(fullfile(p.saveFolder, 'chiasmataSplitting', 'redundancies'), 'redundancies');
    %%
    % NOTE(amotta): Sort queries to make processing easier
    [~, sortedRows] = sortrows(queries(:, 1:3));
    queries = queries(sortedRows, :);
    taskIds = taskIds(sortedRows);
    clear sortedRows;

    % NOTE(amotta): Fix dimensions of `ff` fields
    ff = structfun(@(x) x(:), ff, 'UniformOutput', false);

    axonCount = numel(agglos.axons);
    inputArgs = cell(axonCount, 1);
    sharedInputArgs = {p};

    for curAxonIdx = 1:axonCount
        % find all queries in current axon
        curQueryIds = find(queries(:, 1) == curAxonIdx);
        curQueries = queries(curQueryIds, :);
        assert(issorted(curQueries(:, 1:3), 'rows'));

        % find task IDs and ff indices
        curTaskIds = taskIds(curQueryIds);
        curFfIds = cellfun(@(x) ...
            find(strcmp(ff.filenamesShort, x)), curTaskIds);

        assert(numel(curFfIds) == size(curQueries, 1));
        curChiasmaCount = max(cat(1, 0, curQueries(:, 2)));

        % group flights by chiasma
        curFfGroups = arrayfun( ...
            @(idx) curFfIds(curQueries(:, 2) == idx), ...
            (1:curChiasmaCount)', 'UniformOutput', false);
        curCenterIds = unique(curQueries(:, 7), 'stable');

        % create task defintions
        tasks = struct( ...
            'tracings', cellfun(@(ffIds) struct( ...
                'segids', ff.segIds(ffIds), ...
                'nodes', ff.nodes(ffIds)), ...
                curFfGroups, 'UniformOutput', false), ...
            'centeridx', num2cell(curCenterIds));

        agglos.axons(curAxonIdx).nodesScaled = bsxfun(@times, ...
            agglos.axons(curAxonIdx).nodes(:, 1:3), p.raw.voxelSize);
        inputArgs{curAxonIdx} = {agglos.axons(curAxonIdx), tasks,redundancies{curAxonIdx}};
    end

    job = Cluster.startJob( ...
        @connectEM.splitChiasmataMulti, inputArgs, ...
        'name', mfilename(), ...
        'sharedInputs', sharedInputArgs, ...
        'taskGroupSize', ceil(axonCount / 500), ...
        'numOutputs', 2);
    Cluster.waitForJob(job);
    
    %% save results
    data = fetchOutputs(job);
    data = cat(1, data{:});
    
    % NOTE(amotta): Add `endings` field to old agglomerates
   [oldAgglos.axons.endings] = deal([]);
    oldAgglos.indBigAxons = oldAgglos.indBigAxons(:);
    
    out = struct;
    out.p = p;
    out.oldAgglos = oldAgglos;
    out.gitInfo = Util.gitInfo();
    
    out.summary = cat(1, data{:, 2});
    out.summaryIds = find(oldAgglos.indBigAxons);
    
    out.newAgglos = cat(1, data{:, 1});
    out.parentIds = repelem( ...
        out.summaryIds, cellfun(@numel, data(:, 1)));
    
    % add small agglomerates
    out.newAgglos = cat( ...
        1, out.newAgglos, oldAgglos.axons(~oldAgglos.indBigAxons));
    out.parentIds = cat( ...
        1, out.parentIds, find(~oldAgglos.indBigAxons));
    
    outFile = sprintf('%s-results.mat', datestr(now, 30));
    outFile = fullfile(outputDir, outFile);
    Util.saveStruct(outFile, out);
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