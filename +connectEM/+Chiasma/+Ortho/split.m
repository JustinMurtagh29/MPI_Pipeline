function out = split(param, chiParam, oldAxons, tasks, taskIds, queries)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % restrict to large axons
    axonIds = find(oldAxons.indBigAxons(:));
    axons = oldAxons.axons(axonIds);

    %% replace `nmlFile` by `taskId`
    [~, tasks.id] = ismember(tasks.nmlFile, taskIds.nmlFile);
    tasks.id = tasksIds.id(tasks.id);

    tasks.nmlFile = [];
    tasks = circshift(tasks, 1, 2);
    
    %% add task data to queries
    [~, taskRow] = ismember(queries.taskId, tasks.id);
    queries = cat(2, queries, tasks(taskRow, :));
    queries(:, {'id', 'taskId'}) = [];
    clear taskRow;

    %% statistics
    treeMask = queries.nrInvalidTrees > 0;
    commentMask = queries.nrInvalidComments > 0;
    exitMismatchMask = cellfun( ...
        @(a, b) size(a, 1) ~= size(b, 1), ...
        queries.exits, queries.exitNodeIds);

    fprintf('\n');
    fprintf('# queries answered: %d\n', size(queries, 1));
    fprintf('# queries with invalid tree: %d\n', sum(treeMask));
    fprintf('# queries with invalid comment: %d\n', sum(commentMask));
    fprintf('# queries with invalid exits: %d\n', sum(exitMismatchMask));

    fprintf('⇒ Removing these queries...\n');
    queries(treeMask | commentMask | exitMismatchMask, :) = [];
    clear commentMask missingExitMask;

    fprintf('\n');
    fprintf('# queries left: %d\n', size(queries, 1));

    % sort by axon
    queries = sortrows(queries, 'axonId');

    %% apply results
    uniAxonIds = unique(queries.axonId);
    uniAxonCount = numel(uniAxonIds);
    axonsSplit = cell(uniAxonCount, 1);

    for curIdx = 1:uniAxonCount
        curAxonId = uniAxonIds(curIdx);
        curMask = queries.axonId == curAxonId;

        axonsSplit{curIdx} = ...
            connectEM.Chiasma.Ortho.splitWithQueries( ...
                param, chiParam, axons(curAxonId), queries(curMask, :));
    end

    %% build output
    out = struct;
    out.info = info;

    out.axons = cat(1, axonsSplit{:});
    out.parentIds = repelem( ...
        axonIds(uniAxonIds), cellfun(@numel, axonsSplit));
    otherAxons = setdiff((1:numel(oldAxons.axons))', out.parentIds);

    % add small agglomerates
    out.axons = cat(1, out.axons, oldAxons.axons(otherAxons));
    out.parentIds = cat(1, out.parentIds, otherAxons);

    % build `indBigAxons` mask
    out.indBigAxons = oldAxons.indBigAxons(out.parentIds);
end