function nucleiNew = attachNuclei(graph, nuclei, dendrites)
    % Resort nuclei to match dendrites

    % Find pairwise proabability between each dendrite and axon component
    nucleiNeighbours = cellfun(@(x)unique(cat(2, graph.neighbours{x})), nuclei, 'uni', 0);
    numberNeighbours = zeros(numel(nuclei), numel(dendrites));
    for i=1:numel(dendrites)
        numberNeighbours(:,i) = cellfun(@(x)sum(ismember(x, dendrites{i})), nucleiNeighbours);
    end

    % Find maximal maximum number of neighbours to each dendrite
    [maxNeighbours, maxIdx] = max(numberNeighbours, [], 2);
    isAttachted = maxNeighbours > 10 & maxIdx < 150;

    nucleiNew = cell(max(maxIdx(isAttachted)), 1);
    nucleiNew(maxIdx(isAttachted)) = nuclei(isAttachted);

end

