function splitChiasmataMultiSuper(p)
    % Written by
    %   Kevin Boergens <kevin.boergens@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    oldAxons = load(fullfile( ...
        p.saveFolder, 'aggloState', 'axons_06_c.mat'));
    axons = oldAxons.axons(oldAxons.indBigAxons);
    
    curDir = fullfile( ...
        p.saveFolder, 'chiasmataSplitting', ...
        '20171006T113613-kmb-on-axons-6c');
    
    % set and create output directory
    outputDir = fullfile(curDir, 'outputs');
    if ~exist(outputDir, 'dir'); mkdir(outputDir); end

    load(fullfile(curDir, 'input-data.mat'));
    clear curDir;

    %%
    % NOTE(amotta): Sort queries to make processing easier
    [~, sortedRows] = sortrows(queries(:, 1:3));
    queries = queries(sortedRows, :);
    taskIds = taskIds(sortedRows);
    clear sortedRows;

    % NOTE(amotta): Fix dimensions of `ff` fields
    ff = structfun(@(x) x(:), ff, 'UniformOutput', false);

    axonCount = numel(axons);
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

        axons(curAxonIdx).nodesScaled = bsxfun(@times, ...
            axons(curAxonIdx).nodes(:, 1:3), p.raw.voxelSize);
        inputArgs{curAxonIdx} = {axons(curAxonIdx), tasks};
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
   [oldAxons.axons.endings] = deal([]);
    oldAxons.indBigAxons = oldAxons.indBigAxons(:);
    
    out = struct;
    out.p = p;
    out.oldAxons = oldAxons;
    out.gitInfo = Util.gitInfo();
    
    out.summary = cat(1, data{:, 2});
    out.summaryIds = find(oldAxons.indBigAxons);
    
    out.axons = cat(1, data{:, 1});
    out.parentIds = repelem( ...
        out.summaryIds, cellfun(@numel, data(:, 1)));
    
    % add small agglomerates
    out.axons = cat( ...
        1, out.axons, oldAxons.axons(~oldAxons.indBigAxons));
    out.parentIds = cat( ...
        1, out.parentIds, find(~oldAxons.indBigAxons));
    
    % build `indBigAxons` mask
    out.indBigAxons = oldAxons.indBigAxons(out.parentIds);
    
    outFile = sprintf('%s-results.mat', datestr(now, 30));
    outFile = fullfile(outputDir, outFile);
    Util.saveStruct(outFile, out);
end