function summaries = runFullyQueried()
    % This is a modified version of `connectEM.splitChiasmataMultiSuper`
    % from commit c77836bd23ed5fa37144384da669fc4c076f9cd8.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>

    %% configuration
    rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
    axonFile = fullfile(rootDir, 'aggloState', 'axons_06_c.mat');

    curDir = fullfile( ...
        rootDir, 'chiasmataSplitting', ...
        '20171009T193744-kmb-on-axons-6c');

    % List with input files **in decreasing order of dominance**.
    % Add new requery rounds to the top of this list.
    dataFiles = { ...
        'requeries/20171027T102403_input-data.mat';
        'requeries/20171023T102000_input-data.mat';
        '20171018T205736_input-data.mat'};
    dataFiles = fullfile(curDir, dataFiles);
    clear curDir;

    %% load data
    load(fullfile(rootDir, 'allParameter.mat'), 'p');
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
    axons = oldAxons.axons(oldAxons.indBigAxons);

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

    % Limit ourselves to done chiasmata
    queries = queries(ismember( ...
        queries.uniChiasmaId, uniChiasmaDoneIds), :);

    queries.flightNodes = flights.nodes(queries.flightId);
    queries.flightSegIds = flights.segIds(queries.flightId);
    queries.flightComment = flights.comments(queries.flightId);

    % removing invalid flights and their chiasmata
    invalidFlightMask = ...
        cellfun(@isempty, queries.flightNodes) ...
     | ~cellfun(@isempty, queries.flightComment);

    invalidChiasmaIds = unique(queries.uniChiasmaId(invalidFlightMask));
    queries(ismember(queries.uniChiasmaId, invalidChiasmaIds), :) = [];

    [~, ~, queries.uniChiasmaId] = unique( ...
        queries(:, {'axonId', 'chiasmaId'}), 'rows');

    % Clean-up
    queries(:, {'flightId', 'uniChiasmaId'}) = [];

    %% process axons
    % Group queries by axons
    queries.execute = false(size(queries, 1), 1);
    queries.row = reshape(1:size(queries, 1), [], 1);
    [uniAxonIds, ~, queries.uniAxonId] = unique(queries.axonId);

    % Start applying chiasmata
    queries.overlaps(:) = {[]};
    summaries = cell(numel(uniAxonIds), 1);

    tic
    for curAxonIdx = 1:numel(uniAxonIds)
        curAxon = axons(uniAxonIds(curAxonIdx));
        curAxon.nodesScaled = bsxfun(@times, ...
            curAxon.nodes(:, 1:3), p.raw.voxelSize);
        curQueries = queries( ...
            queries.uniAxonId == curAxonIdx, :);

       [~, curSummary] = ...
            connectEM.splitChiasmataMulti( ...
                p, curAxon, curQueries, 'dryRun', true);

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
end