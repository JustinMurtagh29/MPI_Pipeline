function skel = evaluateAggloCleanSkel(skel, idx, criterium)
    while true
        adM = skel.createAdjacencyMatrix(idx);
        candidates = find((sum(adM, 2) == 1) & criterium);
        if isempty(candidates)
            break;
        end

        skel = skel.deleteNodes(idx, candidates, false);
        criterium(candidates) = [];
        assert(graphconncomp(sparse(skel.createAdjacencyMatrix(idx))) == 1);
    end
end
