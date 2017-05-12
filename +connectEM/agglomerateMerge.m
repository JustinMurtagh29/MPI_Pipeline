function [axonsNew, unchangedIdx] = agglomerateMerge(graph, segmentMeta, borderMeta, axons, result);

    options.latentScore = 0.8;
    options.segDirScore = 0.9;
    options.neuriCScore = 0.7;

    % Create lookup of agglo index for all segments
    axonsLookup = createLookup(segmentMeta, axons);
    % Merge all endings that touch each other (given some criteria)
    idxDirectional = cellfun(@(x)x > options.latentScore, result.latent);
    idxEnding = cellfun(@(x)abs(x) > options.segDirScore, result.scores(idxDirectional), 'uni', 0);
    endingEdges = cellfun(@(x,y)x(y,:), results.edges(idxDirectional), idxEnding, 'uni', 0); 
    endingProb = cellfun(@(x,y)x(y), results.prob(idxDirectional), idxEnding, 'uni', 0); 
    endingBorderIdx = cellfun(@(x,y)x(y), results.borderIdx(idxDirectional), idxEnding, 'uni', 0); 

end

function agglos_reverse = createLookup(segmentMeta, agglos)

    agglos_reverse = zeros(size(segmentMeta.point, 1), 1);
    for idx = 1 : length(agglos)
        agglos_reverse(agglos{idx}) = idx;
    end

end

