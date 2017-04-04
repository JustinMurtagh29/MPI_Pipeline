function [dendrites, axons, er] = extractAndTransferER(graph, dendrites, axons)

    % Find pairwise proabability between each dendrite and axon component
    axonNeighbours = cellfun(@(x)unique(cat(2, graph.neighbours{x})), axons, 'uni', 0);
    axonNeighProb = cellfun(@(x)unique(cat(1, graph.neighProb{x})), axons, 'uni', 0);
    probabilities = zeros(numel(dendrites), numel(axons));
    for i=1:numel(dendrites)
        idx = cellfun(@(x)ismember(x, dendrites{i}), axonNeighbours, 'uni', 0);
        probabilities(i,:) = cellfun(@(x,y)sum(x(y)), axonNeighProb, idx);
    end

    % Find maximal probability of each axon, reassign if above 200% summed
    [maxProb, maxIdx] = max(probabilities, [], 1);
    isEr = maxProb > 2;
    er = axons(isEr);
    assignTo = maxIdx(isEr);

    % Remove ER components from axon and assign to respective dendrite class
    axons = axons(~isEr);
    for i=1:length(er)
        dendrites{assignTo(i)} = cat(1, dendrites{assignTo(i)}, er{i});
    end

end

