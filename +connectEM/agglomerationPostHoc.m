function todo = agglomerationPostHoc(gridAgglo_05, graph, excludedSegmentIdx, segmentMeta, agnostic)
    function donothing()
    end
    todo = [];
    for idx = 1 : length(gridAgglo_05{564}.axonsFinal)
        if ~mod(idx, 100)
            idx
        end
        thisagglo = gridAgglo_05{564}.axonsFinal{idx};
        allNeighbours= graph.neighbours(thisagglo);
        allNeighboursFlat = cell2mat(allNeighbours')';

        outnodes = setdiff(allNeighboursFlat, thisagglo);
        outnodes(excludedSegmentIdx(outnodes)) = [];
        if ~agnostic
            outnodes(segmentMeta.axonProb(outnodes)<0.5) = [];
        end
        allNeighProbs = graph.neighProb(thisagglo);
        allNeighProbsFlat = cell2mat(allNeighProbs);
        neighbourCount = cellfun(@length, allNeighbours);
        allNeighIdxFlat = repelem(1 : length(thisagglo), neighbourCount);
        allNeighFlat = [allNeighboursFlat, allNeighProbsFlat, allNeighIdxFlat'];
        neighPerOutnode = ...
            allNeighFlat(ismember(allNeighboursFlat, outnodes), :);
        [~, maxProbOutnodeIdx] = max(neighPerOutnode(:, 2));
        if isempty(maxProbOutnodeIdx)
            continue
        end
        segmentIdxWithinAgglo = ...
            thisagglo(neighPerOutnode(maxProbOutnodeIdx, 3));
        segmentIdxOfOutnode = ...
            neighPerOutnode(maxProbOutnodeIdx, 1);
        todo(end + 1, :) = sort([segmentIdxWithinAgglo, segmentIdxOfOutnode]);
        segmentMeta.voxelCount(segmentIdxOfOutnode) = Inf;
        % for idx3 = 1 : length(outnodes)
        %     tannstep1 = cellfun(@max, graph.neighProb(thisagglo));
        %
        %     tannpre{idx3} = cellfun(@(x, y) max([-1; max(y(x==outnodes(idx3)))]), graph.neighbours(thisagglo),graph.neighProb(thisagglo), 'uni', 0);
        %     tann(idx3) = max(cell2mat(tannpre{idx3}));
        % end
        % [~, I] = max(tann);
        % [~, I2] = max(cell2mat(tannpre{I}));
        %
        % todo(end + 1, :) = sort(outnodes(I), thisagglo(I2));
    end
end
