function out = splitChiasmataMultiSuper(p, axonFile, dataFiles)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>

    %% load agglomerates
    oldAxons = load(axonFile);

    % set default value for `endings`
    if ~isfield(oldAxons.axons, 'endings')
       [oldAxons.axons.endings] = deal([]);
    end

    % set default value for `solvedChiasma`
    if ~isfield(oldAxons.axons, 'solvedChiasma')
        solvedChiasma = arrayfun( ...
            @(a) false(size(a.nodes, 1), 1), ...
            oldAxons.axons, 'UniformOutput', false);
       [oldAxons.axons.solvedChiasma] = deal(solvedChiasma{:});
        clear solvedChiasma;
    end

    oldAxons.indBigAxons = oldAxons.indBigAxons(:);
    bigAxonIds = find(oldAxons.indBigAxons);
    axons = oldAxons.axons(bigAxonIds);

    %% load queries / flights
    dataFileCount = numel(dataFiles);
    queries = cell(dataFileCount, 1);
    flights = cell(dataFileCount, 1);

    for curFileIdx = 1:numel(dataFiles)
        curDataFile = dataFiles{curFileIdx};
        curData = load(curDataFile);

        % original queries
        curQueries = table;
        curQueries.axonId = curData.queries(:, 1);
        curQueries.chiasmaId = curData.queries(:, 2);
        curQueries.exitId = curData.queries(:, 3);
        curQueries.exitNodeId = curData.queries(:, 8);
        curQueries.seedPos = curData.queries(:, 4:6);
        curQueries.centerNodeId = curData.queries(:, 7);
        curQueries.taskId = curData.taskIds;
        curFlights = curData.ff;

        queries{curFileIdx} = curQueries;
        flights{curFileIdx} = curFlights;
        clear curData;
    end

    queries = cat(1, queries{:});
    flights = Util.concatStructs('last', flights{:});

    fprintf('# handed out queries: %d\n', size(queries, 1));
    fprintf('# queries answered: %d\n', numel(flights.segIds));

    %% cleaning up input data
    % Remove duplicate queries (e.g., due to requerying)
    % This also sorts the queries to make processing easier.
    [~, uniRows] = unique(queries(:, ...
        {'axonId', 'chiasmaId', 'exitId'}), 'rows');

    queries = queries(uniRows, :);
    flights = structfun(@(x) x(:), flights, 'UniformOutput', false);

    % Assign flight paths to queries
    [~, queries.flightId] = ismember( ...
        queries.taskId, flights.filenamesShort);

    % Find chiasmata for which we have all queries answered
    [~, ~, queries.uniChiasmaId] = unique( ...
        queries(:, {'axonId', 'chiasmaId'}), 'rows');
    uniChiasmaDoneIds = find(accumarray( ...
        queries.uniChiasmaId, queries.flightId, [], @all));
    fprintf('# chiasmata answered: %d\n', numel(uniChiasmaDoneIds));

    % Limit ourselves to done chiasmata
    queries = queries(ismember( ...
        queries.uniChiasmaId, uniChiasmaDoneIds), :);

    queries.flightNodes = flights.nodes(queries.flightId);
    queries.flightSegIds = flights.segIds(queries.flightId);
    queries.flightComment = flights.comments(queries.flightId);

    fprintf('\n');
    fprintf('# queries without nodes: %d\n', ...
        sum(cellfun(@isempty, queries.flightNodes)));
    fprintf('# queries with comment: %d\n', ...
        sum(~cellfun(@isempty, queries.flightComment)));
    fprintf('⇒ Removing the above queries ...\n');

    % removing invalid flights and their chiasmata
    invalidFlightMask = ...
        cellfun(@isempty, queries.flightNodes) ...
     | ~cellfun(@isempty, queries.flightComment);

    invalidChiasmaIds = unique(queries.uniChiasmaId(invalidFlightMask));
    queries(ismember(queries.uniChiasmaId, invalidChiasmaIds), :) = [];

    [~, ~, queries.uniChiasmaId] = unique( ...
        queries(:, {'axonId', 'chiasmaId'}), 'rows');
    uniChiasmaDoneIds = find(accumarray( ...
        queries.uniChiasmaId, queries.flightId, [], @all));

    fprintf('\n');
    fprintf('# answered chiasmata left: %d\n', numel(uniChiasmaDoneIds));

    % Clean-up
    queries(:, {'flightId', 'uniChiasmaId'}) = [];

    %% process axons
    % Group queries by axons
    queries.execute = false(size(queries, 1), 1);
    queries.row = reshape(1:size(queries, 1), [], 1);
    [uniAxonIds, ~, queries.uniAxonId] = unique(queries.axonId);

    % Start applying chiasmata
    queries.overlaps(:) = {[]};

    axonsSplit = cell(numel(uniAxonIds), 1);
    summaries = cell(numel(uniAxonIds), 1);

    tic
    for curAxonIdx = 1:numel(uniAxonIds)
        curAxon = axons(uniAxonIds(curAxonIdx));
        curAxon.nodesScaled = bsxfun(@times, ...
            curAxon.nodes(:, 1:3), p.raw.voxelSize);
        curQueries = queries( ...
            queries.uniAxonId == curAxonIdx, :);

       [axonsSplit{curAxonIdx}, curSummary] = ...
            connectEM.splitChiasmataMulti(p, curAxon, curQueries);

        curOverlaps = cellfun(@(t) ...
            t.overlaps, curSummary.tracings, ...
            'UniformOutput', false);
        curOverlaps = cat(1, curOverlaps{:});

        summaries{curAxonIdx} = curSummary;
        queries.overlaps(curQueries.row) = curOverlaps;
    end
    toc

    summaries = cat(1, summaries{:});
    queries.uniAxonId = [];

    %% evaluate flight path attachment
    numEndings = cellfun(@(o) sum(o > 0), queries.overlaps);

    flightEval = table;
    flightEval.twoEndings = (numEndings == 2);
    flightEval.tooFewEndings = (numEndings < 2);

    temp = varfun(@sum, flightEval);
    temp.Properties.VariableNames = flightEval.Properties.VariableNames;

    fprintf('\n');
    fprintf('Flight path evaluation:\n\n');
    disp(temp);

    %% evaluate partitions
    chiasma = arrayfun(@(s) ...
        cellfun(@(t) { ...
            t.partition(:), t.partitionIsValid}, ...
            s.tracings, 'UniformOutput', false), ...
        summaries, 'UniformOutput', false);
    chiasma = cat(1, chiasma{:});
    chiasma = cat(1, chiasma{:});

    chiasma = cell2table( ...
        chiasma, 'VariableNames', {'partition', 'valid'});

    makePartitionStr = @(p) strjoin( ...
        arrayfun(@num2str, p(:)', 'Uni', false), '-');
    chiasmaPartitionStr = cellfun( ...
        makePartitionStr, chiasma.partition, 'UniformOutput', false);
    [uniPartitions, ~, uniPartitionIds] = unique(chiasmaPartitionStr);

    partitionEval = table;
    partitionEval.partition = uniPartitions;
    partitionEval.total = accumarray(uniPartitionIds, 1);
    partitionEval.valid = accumarray(uniPartitionIds, chiasma.valid);
    partitionEval.invalid = partitionEval.total - partitionEval.valid;
    partitionEval = partitionEval(:, ...
        {'partition', 'valid', 'invalid', 'total'});

    fprintf('Evaluation of chiasmata partitioning:\n\n');
    disp(partitionEval);

    %% evaluate solved chiasmata
    chiEval = table;
    chiEval.split = cat(1, summaries.split);
    chiEval.solved = cat(1, summaries.solved);

    fprintf('# chiasmata answered: %d\n', size(chiEval, 1));
    fprintf('# chiasmata marked for splitting: %d\n', sum(chiEval.split));
    fprintf('# chiasmata actually solved: %d\n', sum(chiEval.solved));

    %%
    out = struct;
    out.info = info;

    out.p = p;
    out.oldAxons = oldAxons;

    out.summary = summaries;
    out.summaryIds = bigAxonIds(uniAxonIds);

    out.axons = cat(1, axonsSplit{:});
    out.parentIds = repelem(out.summaryIds, cellfun(@numel, axonsSplit));
    otherAxons = setdiff((1:numel(oldAxons.axons))', out.summaryIds);

    % add small agglomerates
    out.axons = cat(1, out.axons, oldAxons.axons(otherAxons));
    out.parentIds = cat(1, out.parentIds, otherAxons);

    % build `indBigAxons` mask
    out.indBigAxons = oldAxons.indBigAxons(out.parentIds);
end
