function visualizeSuperaggloEvolution(idxState1, idxState2, outputFolder)
    % Compares axon superagglos between 2 different states including some statistics and visualization

    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    % Load parameter struct
    load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat', 'p');
    load(fullfile(p.saveFolder, 'segmentMeta.mat'), 'maxSegId');
    state1 = load(fullfile(p.saveFolder, 'aggloState', ['axons_' num2str(idxState1, '%.2i') '.mat']));
    % The reshape necessary as superagglos currently change dimensionality going from step 3->4
    state1 = reshape(state1.axons,[],1);
    state2 = load(fullfile(p.saveFolder, 'aggloState', ['axons_' num2str(idxState2, '%.2i') '.mat']));
    state2 = reshape(state2.axons,[],1);

    % Display statistics on single superagglos
    display(['-- Statistics superagglo state ' num2str(idxState1)]);
    displaySuperaggloStats(state1);
    display(['-- Statistics superagglo state ' num2str(idxState2)]);
    displaySuperaggloStats(state2);

    % Build lookup for each agglo
    lookup1 = buildLookup(maxSegId, state1);
    lookup2 = buildLookup(maxSegId, state2);
    lookup = cat(2, lookup1, lookup2);
    % Exclude added and removed segments for merger and split calculations
    lookupPersistent = lookup;
    lookupPersistent(any(lookupPersistent == 0, 2),:) = [];
    % Keep only unique rows (used to remove redundant segment overlaps)
    lookupPersistent = unique(lookupPersistent, 'rows');


    % Display statistics on difference between states
    display(['-- Statistics superagglo difference between state ' num2str(idxState1) ' & ' num2str(idxState2)]);
    [agglosMerged, agglosSplit, addedSegId, removedSegId] = displaySuperaggloDiffStats(lookup, lookupPersistent);

    % Write .nmls of differences to inspect in wK
    rng default; % Make sure seed is the same every time
    visualizeSplitsAndMergerAsNml(state1, state2, lookupPersistent, agglosMerged, agglosSplit, outputFolder);
    visualizeAddedAndRemovedSegmentsAsNml(state1, state2, lookup, addedSegId, removedSegId, outputFolder);

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
    % Merger show by double occurences in second column of edges
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
    theseAgglos = agglosMerged(idx);
    for i=1:length(theseAgglos)
        agglosBefore = lookupPersistent(lookupPersistent(:,2) == theseAgglos(i),1);
        visualizeSetOfSuperagglos(state1, state2, agglosBefore, theseAgglos(i), strcat(outputFolder,'randomMerger',num2str(i,'%.3i'),'.nml'));
    end

    % Visualize largest merger (most agglos before merged into one)
    agglosAfterUnique = unique(lookupPersistent(:,2));
    mergerCount = histc(lookupPersistent(:,2), agglosAfterUnique);
    [~, idx] = max(mergerCount);
    agglosBefore = lookupPersistent(lookupPersistent(:,2) == agglosAfterUnique(idx),1);
    visualizeSetOfSuperagglos(state1, state2, agglosBefore, agglosAfterUnique(idx), strcat(outputFolder,'largestMerger',num2str(i,'%.3i'),'.nml'));

    % Visualize 100 random splits
    idx = randperm(numel(agglosSplit), min(numel(agglosSplit),100));
    theseAgglos = agglosSplit(idx);
    for i=1:length(theseAgglos)
        agglosAfter = lookupPersistent(lookupPersistent(:,1) == theseAgglos(i),2);
        visualizeSetOfSuperagglos(state1, state2, theseAgglos(i), agglosAfter, strcat(outputFolder,'randomSplit',num2str(i,'%.3i'),'.nml'));
    end

    % Visualize largest split (agglo split into most agglos after)
    agglosBeforeUnique = unique(lookupPersistent(:,1));
    splitCount = histc(lookupPersistent(:,1), agglosBeforeUnique);
    [~, idx] = max(splitCount);
    agglosAfter = lookupPersistent(lookupPersistent(:,1) == agglosBeforeUnique(idx),2);
    visualizeSetOfSuperagglos(state1, state2, agglosBeforeUnique(idx), agglosAfter, strcat(outputFolder,'largestSplit',num2str(i,'%.3i'),'.nml'));

end

function visualizeAddedAndRemovedSegmentsAsNml(state1, state2, lookup, addedSegId, removedSegId, outputFolder)

    addedAgglos = lookup((lookup(:,1) == 0) & (lookup(:,2) ~= 0),2);
    addedAgglosUnique = unique(addedAgglos);
    removedAgglos = unique(lookup((lookup(:,2) == 0) & (lookup(:,1) ~= 0),1));
    removedAgglosUnique = unique(removedAgglos);

    % Visualize 100 random agglos with added segments
    idx = randperm(numel(addedAgglosUnique), min(numel(addedAgglosUnique),100));
    for i=1:length(idx)
        treeName = sprintf('after_aggloIdx_%.7i', addedAgglosUnique(idx(i)));
        outputFile = strcat(outputFolder,'addedSegments',num2str(i,'%.3i'),'.nml');
        visualizeSuperaggloWithComments(state2(idx(i)), addedSegId, 'added', treeName, outputFile);
    end

    % Visualize agglo with most added segments
    addedSegmentCount = histc(addedAgglos, addedAgglosUnique);
    [~, idx] = max(addedSegmentCount);
    treeName = sprintf('after_aggloIdx_%.7i', addedAgglosUnique(idx));
    outputFile = strcat(outputFolder,'addedSegmentsMost.nml');
    visualizeSuperaggloWithComments(state2(idx), addedSegId, 'added', treeName, outputFile);

    % Visualize 100 random agglos with removed segments
    idx = randperm(numel(removedAgglosUnique), min(numel(removedAgglosUnique),100));
    for i=1:length(idx)
        treeName = sprintf('before_aggloIdx_%.7i', removedAgglosUnique(idx(i)));
        outputFile = strcat(outputFolder,'removedSegments',num2str(i,'%.3i'),'.nml');
        visualizeSuperaggloWithComments(state1(idx(i)), removedSegId, 'removed', treeName, outputFile);
    end

    % Visualize agglo with most removed segments
    removedSegmentCount = histc(removedAgglos, removedAgglosUnique);
    [~, idx] = max(removedSegmentCount);
    treeName = sprintf('before_aggloIdx_%.7i', removedAgglosUnique(idx));
    outputFile = strcat(outputFolder,'removedSegmentsMost.nml');
    visualizeSuperaggloWithComments(state1(idx), removedSegId, 'removed', treeName, outputFile);

end

function visualizeSetOfSuperagglos(agglosBefore, agglosAfter, beforeIdx, afterIdx, outputFile)

    agglos = cat(1, agglosBefore(beforeIdx), agglosAfter(afterIdx));
    nodes = arrayfun(@(x)x.nodes(:,1:3), agglos, 'uni', 0);
    treeNamesB = arrayfun(@(x)sprintf('before_aggloIdx_%.7i', x), beforeIdx, 'uni', 0);
    treeNamesA = arrayfun(@(x)sprintf('after_aggloIdx_%.7i', x), afterIdx, 'uni', 0);
    treeNames = cat(1, treeNamesB, treeNamesA);
    edges = arrayfun(@(x)x.edges, agglos, 'uni', 0);
    connectEM.generateSkeletonFromNodes(outputFile, nodes, treeNames, [], false, edges);

end

function visualizeSuperaggloWithComments(agglos, segIds, comment, treeNames, outputFile);
% Currently only works for one superagglo

    nodes = arrayfun(@(x)x.nodes(:,1:3), agglos, 'uni', 0);
    theseSegId = arrayfun(@(x)x.nodes(:,4), agglos, 'uni', 0);
    idx = ismember(theseSegId{1}, segIds);
    comments = cell(1);
    comments{1}(idx) = repmat({comment}, sum(idx), 1);
    edges = arrayfun(@(x)x.edges, agglos, 'uni', 0);
    connectEM.generateSkeletonFromNodes(outputFile, nodes, treeNames, comments, false, edges);

end

