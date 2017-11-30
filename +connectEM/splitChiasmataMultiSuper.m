function [out, openExits] = splitChiasmataMultiSuper( ...
        p, chiasmataFile, axonFile, dataFiles, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    opts = struct;
    opts.dryRun = false;
    opts.neuriteType = 'axon';
    opts.partialAnswers = false;

    opts = Util.modifyStruct(opts, varargin{:});
    opts.neuriteType = lower(opts.neuriteType);
    
    %% load all chiasmata
    chiasmata = load(chiasmataFile);
    chiasmaParam = chiasmata.info.param.chiasmaParam;
    chiasmata = chiasmata.chiasmata;
    
    allExits = connectEM.Chiasma.Detect.buildTable(chiasmata);
    allExits = connectEM.Chiasma.Flight.selectExits(allExits, [], inf);
    allExits.Properties.VariableNames{'aggloId'} = 'axonId';
    
    forAllExits = @(f) arrayfun( ...
        @(axonId, chiasmaId, exitId) f(axonId, chiasmaId, exitId), ...
        allExits.axonId, allExits.chiasmaId, allExits.exitId);
    
    allExits.exitNodeId = forAllExits( ...
        @(a, c, e)  chiasmata{a}.queryIdx{c}(e));
    allExits.seedPos = cell2mat(forAllExits( ...
        @(a, c, e) {chiasmata{a}.position{c}(e, :)}));
	allExits.centerNodeId = forAllExits( ...
        @(a, c, ~)  chiasmata{a}.ccCenterIdx(c));
    clear forAllExits;

    %% load queries / flights
    dataFileCount = numel(dataFiles);
    queries = cell(dataFileCount, 1);
    flights = cell(dataFileCount, 1);

    for curFileIdx = 1:numel(dataFiles)
        curDataFile = dataFiles{curFileIdx};
        curData = load(curDataFile);

        % build query table
        curQueries = table;
        curQueries.axonId = curData.queries(:, 1);
        curQueries.chiasmaId = curData.queries(:, 2);
        curQueries.exitId = curData.queries(:, 3);
        curQueries.taskId = curData.taskIds;
        
        curFlights = structfun( ...
            @(field) reshape(field, [], 1), ...
            curData.ff, 'UniformOutput', false);

        queries{curFileIdx} = curQueries;
        flights{curFileIdx} = curFlights;
        clear curData;
    end

    queries = cat(1, queries{:});
    flights = Util.concatStructs(1, flights{:});
    
   [~, chiasmaCount] = unique( ...
       queries(:, {'axonId', 'chiasmaId'}), 'rows');
    fprintf('# chiasmata queried: %d\n', numel(chiasmaCount));
    clear chiasmaCount;
    
    fprintf('# handed out queries: %d\n', size(queries, 1));
    fprintf('# queries answered: %d\n', numel(flights.segIds));

    %% cleaning up the flight paths
    % restrict to valid flights
    fprintf('\n');
    invalidNodeMask = cellfun(@isempty, flights.nodes);
    fprintf('# queries without nodes: %d\n', sum(invalidNodeMask));
    invalidCommentMask = ~cellfun(@isempty, flights.comments);
    fprintf('# queries with comment: %d\n', sum(invalidCommentMask));
    fprintf('⇒ Removing the above queries ...\n');

    % removing invalid flights and their chiasmata
    validFlightMask = ~(invalidNodeMask | invalidCommentMask);
    clear invalidNodeMask invalidCommentMask;
    
    flights = structfun( ...
        @(vals) vals(validFlightMask, :), ...
        flights, 'UniformOutput', false);
    
    %% cleaning up the list of chiasmata
    % Remove duplicate queries (e.g., due to requerying)
    % This also sorts the queries to make processing easier.
   [~, uniRows] = unique(queries(:, ...
        {'axonId', 'chiasmaId', 'exitId'}), 'rows');
    queries = queries(uniRows, :);
    
    % Assign flight paths to queries
   [~, queries.flightId] = ismember( ...
        queries.taskId, flights.filenamesShort);
    
    % Mark exits for which there is no valid answer
    openExits = queries( ...
        ~queries.flightId, {'axonId', 'chiasmaId', 'exitId'});
    openExits.Properties.VariableNames{'axonId'} = 'aggloId';

    % Find chiasmata for which we have all queries answered
   [~, ~, queries.uniChiasmaId] = unique( ...
        queries(:, {'axonId', 'chiasmaId'}), 'rows');
    
    fprintf('\n');
    if opts.partialAnswers
        uniChiasmaIds = find(accumarray( ...
            queries.uniChiasmaId, queries.flightId, [], @any));
        
        fprintf('# chiasmata partially answered: %d\n', numel(uniChiasmaIds));
        fprintf('⇒ Limiting to partially answered chiasmata...\n');
    else
        uniChiasmaIds = find(accumarray( ...
            queries.uniChiasmaId, queries.flightId, [], @all));
        
        fprintf('# chiasmata fully answered: %d\n', numel(uniChiasmaIds));
        fprintf('⇒ Limiting to fully answered chiasmata...\n');
    end
    
    queryMask = ismember(queries.uniChiasmaId, uniChiasmaIds);
    queries = queries(queryMask, :);
    clear uniChiasmaIds queryMask;
    
    % restrict exits
    exitMask = ismember( ...
        allExits(:, {'axonId', 'chiasmaId'}), ...
        queries(:, {'axonId', 'chiasmaId'}), 'rows');
    allExits = allExits(exitMask, :);
    clear exitMask;
    
    %% add queries to exits
   [curMask, curRows] = ismember( ...
       allExits(:, {'axonId', 'chiasmaId', 'exitId'}), ...
       queries(:, {'axonId', 'chiasmaId', 'exitId'}), 'rows');
    curRows = curRows(curMask);
    
    allExits.taskId(curMask) = queries.taskId(curRows);
    allExits.flightId(curMask) = queries.flightId(curRows);
    clear curMask curRows;
    
    %% add flights
    curMask = allExits.flightId > 0;
    curRows = allExits.flightId(curMask);
    
    allExits.flightNodes(curMask) = flights.nodes(curRows);
    allExits.flightSegIds(curMask) = flights.segIds(curRows);
    allExits.flightComment(curMask) = flights.comments(curRows);
    clear curMask curRows;
    
    % TODO(amotta): Improve
    queries = allExits;
    
    %% load agglomerates
    oldAxons = load(axonFile);
    
    switch opts.neuriteType
        case 'axon'
            % nothing to do
        case 'dendrite'
            oldAxons.axons = oldAxons.dendrites;
            oldAxons.indBigAxons = oldAxons.indBigDends;
        otherwise
            error('Unknown neurite type "%s"', opts.neuriteType)
    end
    
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

    %% process axons
    queries.execute(:) = false;
    queries.overlaps(:) = cell(1, 1);

    % Make sure that queries are properly sorted
    % This is required for the assignment of `overlaps` to work.
    assert(issortedrows(queries(:, {'axonId', 'chiasmaId', 'exitId'})));
    
    % Find axons which need to be processed
   [uniAxonIds, ~, queries.uniAxonId] = unique(queries.axonId);

    axonsSplit = cell(numel(uniAxonIds), 1);
    summaries = cell(numel(uniAxonIds), 1);
    
    splitOpts = horzcat(fieldnames(opts), struct2cell(opts));
    splitOpts = reshape(transpose(splitOpts), 1, []);

    tic
    for curAxonIdx = 1:numel(uniAxonIds)
        curAxon = axons(uniAxonIds(curAxonIdx));
        curAxon.nodesScaled = bsxfun(@times, ...
            curAxon.nodes(:, 1:3), p.raw.voxelSize);
        
        curRows = find(queries.uniAxonId == curAxonIdx);
        curQueries = queries(curRows, :);
        
        % split axons and save summary
       [axonsSplit{curAxonIdx}, curSummary] = ...
            connectEM.splitChiasmataMulti( ...
                p, chiasmaParam, curAxon, curQueries, splitOpts{:});
        summaries{curAxonIdx} = curSummary;
        
        % build overlaps
        curOverlaps = cellfun(@(t) ...
            t.overlaps, curSummary.tracings, 'UniformOutput', false);
        curOverlaps = cat(1, curOverlaps{:});
        
        curRows = curRows(curQueries.flightId > 0);
        queries.overlaps(curRows) = curOverlaps;
    end
    toc

    summaries = cat(1, summaries{:});
    queries.uniAxonId = [];

    %% evaluate flight path attachment
    numEndings = cellfun(@(o) sum(o > 0), queries.overlaps);

    flightEval = table;
    flightEval.attached = (numEndings == 2);
    flightEval.dangling = (numEndings == 1);

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

    %% build output
    out = struct;
    out.p = p;
    out.oldAxons = oldAxons;

    out.summary = summaries;
    out.summaryIds = bigAxonIds(uniAxonIds);
    
    if opts.dryRun; return; end
    
    % handle split agglomerates
    out.axons = cat(1, axonsSplit{:});
    out.parentIds = repelem(out.summaryIds, cellfun(@numel, axonsSplit));
    otherAxons = setdiff((1:numel(oldAxons.axons))', out.summaryIds);
    
    % add small agglomerates
    out.axons = cat(1, out.axons, oldAxons.axons(otherAxons));
    out.parentIds = cat(1, out.parentIds, otherAxons);
    
    % build `indBigAxons` mask
    out.indBigAxons = oldAxons.indBigAxons(out.parentIds);
    
    switch opts.neuriteType
        case 'axon'
            % nothing to do
        case 'dendrite'
            out = renamefield(out, 'axons', 'dendrites');
            out = renamefield(out, 'oldAxons', 'oldDendrites');
            out = renamefield(out, 'indBigAxons', 'indBigDends');
        otherwise
            error('Unknown neurite type "%s"', opts.neuriteType)
    end
end

function s = renamefield(s, old, new)
    assert(~isfield(s, new));
    
   [s.(new)] = s.(old);
    s = rmfield(s, old);
end
