function [axonsNew, changedIdx, unchangedResult] = agglomerateMerge(graph, segmentMeta, borderMeta, axons, result, options);

    assert(numel(axons) == numel(result.neighbours));
    % Create lookup of agglo index for all segments
    axonsLookup = createLookup(segmentMeta, axons);

    % Find all endings (given criteria on latent, directionality, continuity and size of border)
    idxDirectional = cellfun(@(x)x(:,1) > options.latentScore, result.latent, 'uni', 0);
    idxEnding = cellfun(@(x)abs(x) > options.segDirScore, result.scores, 'uni', 0);
    idxContinuity = cellfun(@(x)x > options.neuriCScore, result.prob, 'uni', 0);
    idxLarge = cellfun(@(x)borderMeta.borderSize(x) > options.borderSize, result.borderIdx, 'uni', 0);
    idxAll = cellfun(@(w,x,y,z)w&x&y&z, idxDirectional, idxEnding, idxContinuity, idxLarge, 'uni', 0);
    nrEndings = cellfun(@sum, idxAll);
    
    % Source & target agglomerate
    sources = repelem(1:length(nrEndings), nrEndings)';
    targets = cell2mat(cellfun(@(x,y)axonsLookup(x(y)), result.neighbours, idxAll, 'uni', 0));
    % Exclude small sources & 0 targets (not current in axon agglomerates)
    sourceSize = cellfun(@(x)sum(segmentMeta.voxelCount(x)), axons);
    idx = ismember(sources, find(sourceSize < options.sourceSize)) | targets == 0;
    sources(idx) = [];
    targets(idx) = [];

    % Display some statistics of inital detections
    displayStats(idxDirectional, idxEnding, idxContinuity, idxLarge, idxAll, idx);

    % Do the merging
    edgesBetweenAgglos = cat(2, sources, targets);
    cc = Graph.findConnectedComponents(edgesBetweenAgglos);
    axonsNew = cellfun(@(x)cat(1,axons{x}), cc, 'uni', 0);
    unchangedIdx = setdiff(1:length(axons), cat(1, cc{:}));
    axonsNew(numel(cc)+1:numel(cc)+numel(unchangedIdx)) = axons(unchangedIdx);
    unchangedResult = structfun(@(x)x(unchangedIdx), result, 'uni', 0);
    changedIdx = 1:numel(cc)';

end

function agglos_reverse = createLookup(segmentMeta, agglos)
    % Build lookup of agglomerate ID based on segment ID
    agglos_reverse = zeros(size(segmentMeta.point, 1), 1);
    for idx = 1 : length(agglos)
        agglos_reverse(agglos{idx}) = idx;
    end
end

function displayStats(idxDir, idxEnd, idxCon, idxLarge, idxAll, idx)
    % Display some statistics of this agglomeration merging
    display(['# border directional: ' num2str(sum(cellfun(@sum,idxDir)))]);
    display(['# border ending: ' num2str(sum(cellfun(@sum,idxEnd)))]);
    display(['# border continuity: ' num2str(sum(cellfun(@sum,idxCon)))]);
    display(['# border large: ' num2str(sum(cellfun(@sum,idxLarge)))]);
    display(['# endings (all above): ' num2str(sum(cellfun(@sum,idxAll)))]);
    display(['# ending targets not part of agglomerates: ' num2str(sum(idx))]);
end

