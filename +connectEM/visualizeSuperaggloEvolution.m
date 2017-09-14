function visualizeSuperaggloEvolution(idxState1, idxState2, outputFolder,prefix)
    % Compares axon superagglos between 2 different states including some statistics and visualization
    if ~exist('prefix', 'var')
        prefix = 'axons';
    end
    if exist('outputFolder', 'var') && ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    if isnumeric(idxState1)
        fileState1 = [prefix '_' num2str(idxState1, '%.2i') '.mat'];
    else
        fileState1 = idxState1;
    end
    if isnumeric(idxState2)
        fileState2 = [prefix '_' num2str(idxState2, '%.2i') '.mat'];
    else
        fileState2 = idxState2;
    end
    % Load parameter struct
    load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat', 'p');
    load(fullfile(p.saveFolder, 'segmentMeta.mat'), 'maxSegId');
    state1 = load(fullfile(p.saveFolder, 'aggloState', fileState1));
    % The reshape necessary as superagglos currently change dimensionality going from step 3->4
    state1 = reshape(state1.axons,[],1);
    state2 = load(fullfile(p.saveFolder, 'aggloState', fileState2));
    state2 = reshape(state2.axons,[],1);

    % Display statistics on single superagglos
    display(['-- Statistics superagglo state ' fileState1]);
    displaySuperaggloStats(state1);
    display(['-- Statistics superagglo state ' fileState2]);
    displaySuperaggloStats(state2);

    % Build lookup for each agglo
    display('Building lookup table');
    lookup1 = buildLookup(maxSegId, state1);
    lookup2 = buildLookup(maxSegId, state2);
    lookup = cat(2, lookup1, lookup2);
    % Exclude added and removed segments for merger and split calculations
    lookupPersistent = lookup;
    lookupPersistent(any(lookupPersistent == 0, 2),:) = [];
    % Keep only unique rows (used to remove redundant segment overlaps)
    lookupPersistent = unique(lookupPersistent, 'rows');

    % Display statistics on difference between states
    display(['-- Statistics superagglo difference between state ' fileState1 ' & ' fileState2]);
    [agglosMerged, agglosSplit, addedSegId, removedSegId] = displaySuperaggloDiffStats(lookup, lookupPersistent);

    if exist('outputFolder','var')
        % Write .nmls of differences to inspect in wK
        display(['Writing nmls to: ' outputFolder]);
        rng default; % Make sure seed is the same every time
        visualizeSplitsAndMergerAsNml(state1, state2, lookupPersistent, agglosMerged, agglosSplit, outputFolder);
        visualizeAddedAndRemovedSegmentsAsNml(state1, state2, lookup, addedSegId, removedSegId, outputFolder);
    end
end

function displaySuperaggloStats(state)

    segIds = arrayfun(@(x)x.nodes(:,4), state, 'uni', 0);
    segIdsFlat = cat(1, segIds{:});
    segIdsOnly = segIdsFlat(~isnan(segIdsFlat));
    sizeAgglos = sort(cellfun(@numel, segIds), 'descend');

    display(['Number agglomerates: ' num2str(numel(state))]);
    display(['Number segments: ' num2str(numel(segIdsOnly))]);
    display(['Number unique segments: ' num2str(numel(unique(segIdsOnly)))]);
    display(['Number nodes: ' num2str(sum(arrayfun(@(x)size(x.nodes,1), state)))]);
    display(['Number edges: ' num2str(sum(arrayfun(@(x)size(x.edges,1), state)))]);
    sizeString = sprintf('%d ', sizeAgglos(1:10));
    display(['Largest agglomerates (segment-count): ' sizeString]);

end

function [agglosMerged, agglosSplit, addedSegId, removedSegId] = displaySuperaggloDiffStats(lookup, lookupPersistent)

    % Determine merged superagglos
    agglosAfter = unique(lookupPersistent(:,2));
    mergerCount = histc(lookupPersistent(:,2), agglosAfter);
    isMerger = mergerCount > 1;
    agglosMerged = agglosAfter(isMerger);
    display(['Number merged supperagglos: ' num2str(sum(isMerger))]);
    if sum(isMerger)
        display(['Average number of superagglos merged together: ' num2str(mean(mergerCount(isMerger)))]);
        display(['Maximum number of superagglos merged together: ' num2str(max(mergerCount(isMerger)))]);
    end

    % Determine split superagglos
    agglosBefore = unique(lookupPersistent(:,1));
    splitCount = histc(lookupPersistent(:,1), agglosBefore);
    isSplit = splitCount > 1;
    agglosSplit = agglosBefore(isSplit);
    display(['Number split superagglos: ' num2str(sum(isSplit))]);
    if sum(isSplit)
        display(['Average number of a superagglo split into: ' num2str(mean(splitCount(isSplit)))]);
        display(['Maximum number of a superagglo split into: ' num2str(max(splitCount(isSplit)))]);
    end

    % Determine added segments
    addedSegId = find((lookup(:,1) == 0) & (lookup(:,2) ~= 0));
    addedToAgglo = lookup(addedSegId,2);
    addedSegPerAgglo = histc(addedToAgglo, unique(addedToAgglo));
    display(['Added segments: ' num2str(numel(addedSegId))]);
    if ~isempty(addedSegId)
        display(['Average number of segments added per agglomerate: ' num2str(mean(addedSegPerAgglo))]);
        display(['Maximum number of segments added per agglomerate: ' num2str(max(addedSegPerAgglo))]);
    end

    % Determine removed segments
    removedSegId = find((lookup(:,1) ~= 0) & (lookup(:,2) == 0));
    removedFromAgglo = lookup(removedSegId,1);
    removedSegPerAgglo = histc(removedFromAgglo, unique(removedFromAgglo));
    display(['Removed segments: ' num2str(numel(removedSegId))]);
    if ~isempty(removedSegId)
        display(['Average number of segments removed per agglomerate: ' num2str(mean(removedSegPerAgglo))]);
        display(['Maximum number of segments removed per agglomerate: ' num2str(max(removedSegPerAgglo))]);
    end

end

function lookup = buildLookup(maxSegId, superagglos)
    % Analog to buildLUT from Alessandro for superagglos directly
    lookup = zeros(maxSegId, 1);
    agglos = arrayfun(@(x)unique(x.nodes(~isnan(x.nodes(:,4)),4)), superagglos, 'uni', 0);
    lookup(cell2mat(agglos)) = repelem(1:numel(agglos), cellfun(@numel, agglos));
end

function visualizeSplitsAndMergerAsNml(state1, state2, lookupPersistent, agglosMerged, agglosSplit, outputFolder)

    % Visualize 100 random merger
    idx = randperm(numel(agglosMerged), min(numel(agglosMerged), 100));
    agglosAfter = num2cell(agglosMerged(idx));
    agglosBefore = cellfun(@(x)lookupPersistent(lookupPersistent(:,2) == x,1), agglosAfter, 'uni', 0);
    fileNames = arrayfun(@(x)fullfile(outputFolder, ['merger_random_', num2str(x, '%.3i'), '.nml']), 1:numel(idx), 'uni', 0);
    visualizeSetOfSuperagglos(state1, state2, agglosBefore, agglosAfter, fileNames);

    % Visualize largest 10 merger (most agglos before merged into one)
    agglosAfterUnique = unique(lookupPersistent(:,2));
    mergerCount = histc(lookupPersistent(:,2), agglosAfterUnique);
    [~, idx] = sort(mergerCount, 'descend');
    agglosAfter = num2cell(agglosAfterUnique(idx(1:min(sum(mergerCount > 1), 10))));
    agglosBefore = cellfun(@(x)lookupPersistent(lookupPersistent(:,2) == x,1), agglosAfter, 'uni', 0);
    fileNames = arrayfun(@(x)fullfile(outputFolder, ['merger_largest_', num2str(x, '%.3i'), '.nml']), 1:numel(idx), 'uni', 0);
    visualizeSetOfSuperagglos(state1, state2, agglosBefore, agglosAfter, fileNames);

    % Visualize 100 random splits
    idx = randperm(numel(agglosSplit), min(numel(agglosSplit),100));
    agglosBefore = num2cell(agglosSplit(idx));
    agglosAfter = cellfun(@(x)lookupPersistent(lookupPersistent(:,1) == x,2), agglosBefore, 'uni', 0);
    fileNames = arrayfun(@(x)fullfile(outputFolder, ['split_random_', num2str(x, '%.3i'), '.nml']), 1:numel(idx), 'uni', 0);
    visualizeSetOfSuperagglos(state1, state2, agglosBefore, agglosAfter, fileNames);

    % Visualize largest split (agglo split into most agglos after)
    agglosBeforeUnique = unique(lookupPersistent(:,1));
    splitCount = histc(lookupPersistent(:,1), agglosBeforeUnique);
    [~, idx] = sort(splitCount, 'descend');
    agglosBefore = num2cell(agglosBeforeUnique(idx(1:min(sum(splitCount > 1), 10))));
    agglosAfter = cellfun(@(x)lookupPersistent(lookupPersistent(:,1) == x,2), agglosBefore, 'uni', 0);
    fileNames = arrayfun(@(x)strcat(outputFolder, 'split_largest_', num2str(x, '%.3i'), '.nml'), 1:numel(idx), 'uni', 0);
    visualizeSetOfSuperagglos(state1, state2, agglosBefore, agglosAfter, fileNames);

end

function visualizeSetOfSuperagglos(agglosBefore, agglosAfter, beforeIdx, afterIdx, outputFiles)

    for i=1:length(beforeIdx)
        agglos = cat(1, agglosBefore(beforeIdx{i}), agglosAfter(afterIdx{i}));
        nodes = arrayfun(@(x)x.nodes(:,1:3), agglos, 'uni', 0);
        treeNamesB = arrayfun(@(x)sprintf('before_aggloIdx_%.7i', x), beforeIdx{i}, 'uni', 0);
        treeNamesA = arrayfun(@(x)sprintf('after_aggloIdx_%.7i', x), afterIdx{i}, 'uni', 0);
        treeNames = cat(1, treeNamesB, treeNamesA);
        edges = arrayfun(@(x)x.edges, agglos, 'uni', 0);
        connectEM.generateSkeletonFromNodes(outputFiles{i}, nodes, treeNames, [], false, edges);
    end

end

function visualizeAddedAndRemovedSegmentsAsNml(state1, state2, lookup, addedSegId, removedSegId, outputFolder)

    addedAgglos = lookup((lookup(:,1) == 0) & (lookup(:,2) ~= 0),2);
    addedAgglosUnique = unique(addedAgglos);
    removedAgglos = unique(lookup((lookup(:,2) == 0) & (lookup(:,1) ~= 0),1));
    removedAgglosUnique = unique(removedAgglos);

    if ~isempty(addedAgglosUnique)
        % Visualize 100 random agglos with added segments
        idx = randperm(numel(addedAgglosUnique), min(numel(addedAgglosUnique),100));
        treeNames = arrayfun(@(x){sprintf('after_aggloIdx_%.7i', x)}, addedAgglosUnique(idx), 'uni', 0);
        outputFiles = arrayfun(@(x)strcat(outputFolder,'addedSegments_random_',num2str(x,'%.3i'),'.nml'), 1:numel(idx), 'uni', 0);
        visualizeSuperaggloWithComments(state2(idx), addedSegId, 'added', treeNames, outputFiles);

        % Visualize agglo with most added segments
        addedSegmentCount = histc(addedAgglos, addedAgglosUnique);
        [~, idx] = sort(addedSegmentCount, 'descend');
        idx = idx(1:min(numel(idx),10));
        treeNames = arrayfun(@(x){sprintf('after_aggloIdx_%.7i', x)}, addedAgglosUnique(idx), 'uni', 0);
        outputFiles = arrayfun(@(x)strcat(outputFolder,'addedSegments_most_',num2str(x,'%.3i'),'.nml'), 1:numel(idx), 'uni', 0);
        visualizeSuperaggloWithComments(state2(idx), addedSegId, 'added', treeNames, outputFiles);
    end

    if ~isempty(removedAgglosUnique)
        % Visualize 100 random agglos with removed segments
        idx = randperm(numel(removedAgglosUnique), min(numel(removedAgglosUnique),100));
        treeNames = arrayfun(@(x){sprintf('before_aggloIdx_%.7i', x)}, removedAgglosUnique(idx), 'uni', 0);
        outputFiles = arrayfun(@(x)strcat(outputFolder,'removedSegments_random_',num2str(x,'%.3i'),'.nml'), 1:numel(idx), 'uni', 0);
        visualizeSuperaggloWithComments(state1(idx), removedSegId, 'removed', treeNames, outputFiles);

        % Visualize agglo with most removed segments
        removedSegmentCount = histc(removedAgglos, removedAgglosUnique);
        [~, idx] = sort(removedSegmentCount, 'descend');
        idx = idx(1:min(numel(idx),10));
        treeNames = arrayfun(@(x){sprintf('before_aggloIdx_%.7i', x)}, removedAgglosUnique(idx), 'uni', 0);
        outputFiles = arrayfun(@(x)strcat(outputFolder,'removedSegments_most_',num2str(x,'%.3i'),'.nml'), 1:numel(idx), 'uni', 0);
        visualizeSuperaggloWithComments(state1(idx), removedSegId, 'removed', treeNames, outputFiles);
    end

end

function visualizeSuperaggloWithComments(agglos, segIds, comment, treeNames, outputFiles);
    % Currently only works for one superagglo

    for i=1:length(agglos)
        nodes = arrayfun(@(x)x.nodes(:,1:3), agglos(i), 'uni', 0);
        theseSegId = arrayfun(@(x)x.nodes(:,4), agglos(i), 'uni', 0);
        idx = ismember(theseSegId{1}, segIds);
        comments = cell(1);
        comments{1}(idx) = repmat({comment}, sum(idx), 1);
        edges = arrayfun(@(x)x.edges, agglos(i), 'uni', 0);
        connectEM.generateSkeletonFromNodes(outputFiles{i}, nodes, treeNames{i}, comments, false, edges);
    end

end

