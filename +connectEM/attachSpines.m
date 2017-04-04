function [dendritesNew, spinePaths, comment] = attachSpines(graph, segmentMeta, dendrites, axons, maxSteps )
    % Attach spines to a set of dendritic equivalence classes based on shortest path

    % Starting points are spines (not yet collected or excluded) and stop if non dendritic segment reached
    excluded = ~segmentMeta.isDendrite;
    excluded(cat(1, axons{:})) = true;
    targetLabel = zeros(segmentMeta.maxSegId,1);
    for i=1:length(dendrites)
        targetLabel(dendrites{i}) = i;
    end
    spineIds = find(segmentMeta.isSpine & ~excluded);

    % Display some stats
    nrSpineSegments = sum(segmentMeta.isSpine);
    display(['Total #segments classified as spines: ' num2str(nrSpineSegments)]);
    nrSpineSegmentsExcluded = sum(segmentMeta.isSpine & excluded);
    display(['Total #segments not in dendrite class: ' num2str(nrSpineSegmentsExcluded)]);

    % Perform maximum probaility search for each spine head
    spinePaths = cell(numel(spineIds),1);
    comment = cell(numel(spineIds),1);
    for i=1:length(spineIds)
        takenPath = spineIds(i);
        nrSteps = 0;
        neighProb = [];
        neighbours = [];
        while true
            % Update neighbour list based on new agglomerate
            neighProb = cat(1, neighProb, graph.neighProb{takenPath(end)});
            neighbours = cat(2, neighbours, graph.neighbours{takenPath(end)});
            idx = excluded(neighbours);
            neighbours(idx) = [];
            neighProb(idx) = [];
            idx = ~ismember(neighbours, takenPath);
            neighProb = neighProb(idx);
            neighbours = neighbours(idx);
            
            % Add maximum probability segment to path if none of the terminal conditions are met
            if targetLabel(takenPath(end)) ~= 0
                comment{i} = 'attachted';
                break; 
            elseif isempty(neighbours)
                comment{i} = 'no more dendritic neighbours';
                break;
            elseif nrSteps >= maxSteps
                comment{i} = 'maximum steps reached';
                break;
            else
                [maxProb, idx] = max(neighProb);
                takenPath(end+1) = neighbours(idx);
            end
            nrSteps = nrSteps + 1;
        end
        spinePaths{i} = takenPath';
    end

    % Add to dendritic equivalence classes
    idx = find(strcmp(comment, 'attachted'));
    nrSpinesAttachted = numel(idx);
    nrSegmentsCollected = numel(unique(cat(1, spinePaths{idx})));
    display(['Total #segments attachted to dendrite class: ' num2str(nrSpinesAttachted)]);
    display(['Total segments added to dendrite agglomerations in the process: ' num2str(nrSegmentsCollected)]);
    dendritesNew = dendrites;
    for i=1:length(idx)
        dendriteIdx = targetLabel(spinePaths{idx(i)}(end)); 
        dendritesNew{dendriteIdx} = cat(1, dendritesNew{dendriteIdx}, spinePaths{idx(i)});
    end
    dendritesNew = cellfun(@unique, dendritesNew, 'uni', 0);

    % Remove 8 segments that are in more than 1 component
    segIds = cat(1, dendritesNew{:});
    uniqueSegIds = unique(segIds);
    counts = histc(segIds, uniqueSegIds);
    segIdsToExclude = uniqueSegIds(counts > 1);
    dendritesNew = cellfun(@(x)x(~ismember(x, segIdsToExclude)), dendritesNew, 'uni', 0);

end

