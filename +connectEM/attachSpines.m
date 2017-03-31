function [dendritesNew, spinePaths] = attachSpines(graph, segmentMeta, dendrites, maxSteps )
% Attach spines to a set of dendritic equivalence classes based on shortest path

% Starting points are spines (not yet collected or excluded) and stop if non dendritic segment reached
excluded = ~segmentMeta.isDendrite;
targetLabel = zeros(segmentMeta.maxSegId,1);
for i=1:length(dendrites)
    targetLabel(dendrites{i}) = i;
end
spineIds = find(segmentMeta.isSpine & targetLabel == 0 & ~excluded);

% Perform maximum probaility search for each spine head
spinePaths = cell(numel(spineIds),1);
for i=1:length(spineIds)
    takenPath = spineIds(i);
    nrSteps = 1;
    neighProb = [];
    neighbours = [];
    while nrSteps <= maxSteps && targetLabel(takenPath(end)) == 0;
        neighProb = cat(1, neighProb, graph.neighProb{takenPath(end)});
        neighbours = cat(2, neighbours, graph.neighbours{takenPath(end)});
        % Exclude excluded neighbours
        idx = excluded(neighbours);
        neighbours(idx) = [];
        neighProb(idx) = [];
        % Exclude already visited
        idx = ~ismember(neighbours, takenPath);
        neighProb = neighProb(idx);
        neighbours = neighbours(idx);
        % Add maximum probability segment to path if there is still a neighbour around
        [maxProb, idx] = max(neighProb);
        if ~isempty(neighbours) && maxProb > 0.5
            takenPath(end+1) = neighbours(idx);
        else
            break;
        end
        nrSteps = nrSteps + 1;
    end
    if targetLabel(takenPath(end)) ~= 0
        spinePaths{i} = takenPath';
    end
end

% Add to dendritic equivalence classes
idx = find(~cellfun(@isempty, spinePaths));
dendritesNew = dendrites;
for i=1:length(idx)
    dendriteIdx = targetLabel(spinePaths{idx(i)}(end)); 
    dendritesNew{dendriteIdx} = cat(1, dendritesNew{dendriteIdx}, spinePaths{idx(i)});
end
dendritesNew = cellfun(@unique, dendritesNew, 'uni', 0);

end

