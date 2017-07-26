function [axonsNew, dendritesNew, garbageEdges] = garbageCollection(graph, segmentMeta, axons, dendrites, heuristics)

    % Give axons and dendrites same treatment here
    nrAxonComponents = numel(axons);
    theseComponents = cat(1, axons, dendrites);

    % Keep only high probability edges that are not in heuristics
    idx = graph.prob > 0.99 & ~any(ismember(graph.edges, cat(1, heuristics{:})), 2);
    edges = graph.edges(idx,:);
    prob = graph.prob(idx);
    % Keep only the ones with exactly one collectedSegment and have collected segment in first row
    collectedSegments = cat(1, theseComponents{:});
    edgeIdx = ismember(edges, collectedSegments);
    edgeIdx(all(edgeIdx, 2),:) = false;
    [row, col] = find(edgeIdx);
    edges = cat(2, edges(sub2ind(size(edges), row, col)), edges(sub2ind(size(edges), row, abs(col-3))));
    prob = prob(row);
    % Construct target labels for faster execution
    targetLabel = zeros(segmentMeta.maxSegId,1);
    for i=1:length(theseComponents)
        targetLabel(theseComponents{i}) = i;
    end
    % Sort according to probability
    [prob, idx] = sort(prob, 'descend');
    edges = edges(idx,:);
    clear row col idx i;

    initialSize = size(edges,1);
    alreadyCollected = [];
    garbageEdges = zeros(0,2);
    for i=1:initialSize
        if ~ismember(edges(i,2), alreadyCollected)
            theseComponents{targetLabel(edges(i,1))}(end+1) = edges(i,2);
            alreadyCollected(end+1) = edges(i,2);
            garbageEdges(end+1,:) = edges(i,:);
        end
    end

    axonsNew = theseComponents(1:nrAxonComponents);
    dendritesNew = theseComponents(nrAxonComponents+1:end);

end


