function [axonsNew, unchangedIdx] = agglomerateMerge(graph, segmentMeta, borderMeta, axons, result);

    % Thresholds for ending detection
    options.latentScore = 0.8;
    options.segDirScore = 0.9;
    options.neuriCScore = 0.7;
    options.borderSize = 40;

    % Thresholds on each ending pair
    options.maxAngle = 45;
    options.maxSize = 1e7;
    options.probSurplus = 1;

    % Create lookup of agglo index for all segments
    axonsLookup = createLookup(segmentMeta, axons);
    % Find all endings (given criteria on latent, directionality, continuity and size of border)
    idxDirectional = cellfun(@(x)x(:,1) > options.latentScore, result.latent, 'uni', 0);
    idxEnding = cellfun(@(x)abs(x) > options.segDirScore, result.scores, 'uni', 0);
    idxContinuity = cellfun(@(x)x > options.neuriCScore, result.prob, 'uni', 0);
    idxLarge = cellfun(@(x)borderMeta.borderSize(x) > options.borderSize, result.borderIdx, 'uni', 0);
    idxAll = cellfun(@(w,x,y,z)w&x&y&z, idxDirectional, idxEnding, idxContinuity, idxLarge, 'uni', 0);
    displayStats(idxDirectional, idxEnding, idxContinuity, idxLarge, idxAll);
    nrEndings = cellfun(@numel, idxAll);
    % Posititon and direction of ending based on local surround and source & target agglomerate
    sources = repelem(1:length(nrEndings), nrEndings);
    targets = cell2mat(cellfun(@(x,y)axonsLookup(x(y)), result.neighbours, idxAll, 'uni', 0));
    edgesBetweenAgglos = cat(1, sources, targets);
    position = cell2mat(cellfun(@(x,y)borderMeta.borderCoM(x(y,:)), result.borderIdx, idxAll, 'uni', 0));
    direction = cell2mat(cellfun(@(x,y)squeeze(x(:,1,y))', result.pca, idxAll, 'uni', 0));
    % Only merge edges that are detected as endings on both sides
    uniqueEdges = unique(edgesBetweenAgglos, 'rows');
    % Or if one is very small

    % Check angles between local directions

    % Check total size after merging

    % Check probabilities between agglomerates

    % Display some statistics
    displayStats(idxDirectional, idxEnding, idxContinuity, idxLarge, idxAll);
    % Do the merging
    cc = Graph.findConnectedComponents(edgesBetweenAgglos);
    axonsNew = axons{cc};

end

function agglos_reverse = createLookup(segmentMeta, agglos)
    % Build lookup of agglomerate ID based on segment ID
    agglos_reverse = zeros(size(segmentMeta.point, 1), 1);
    for idx = 1 : length(agglos)
        agglos_reverse(agglos{idx}) = idx;
    end
end

function displayStats(idxDir, idxEnd, idxCon, idxLarge, idxAll)
    % Display some statistics of this agglomeration merging
    display(['# border directional: ' num2str(sum(cellfun(@sum,idxDir)))]);
    display(['# border ending: ' num2str(sum(cellfun(@sum,idxEnd)))]);
    display(['# border continuity: ' num2str(sum(cellfun(@sum,idxCon)))]);
    display(['# border large: ' num2str(sum(cellfun(@sum,idxLarge)))]);
    display(['# endings (all above): ' num2str(sum(cellfun(@sum,idxAll)))]);
end

