function gtConsolidated = consolidateRedundantAnnotations( gt )
% Consolidate redundant annotations as extracted by
% getContinuityLabelsFromNml

gtConsolidated = struct();
for i=1:size(gt,1)
    allEdges = cat(1,gt(i,:).edges);
    allLabels = cat(1,gt(i,:).labels);
    allProb = cat(1,gt(i,:).prob);
    [uniqueEdges, ~, idx] = unique(allEdges, 'rows');
    gtConsolidated(i).edges = uniqueEdges;
    for j=1:size(uniqueEdges,1)
        theseIdx = idx == j;
        labels = allLabels(theseIdx);
        if ~isempty(labels(labels ~= 0)) && all(labels(labels ~= 0) == 1)
            gtConsolidated(i).labels(j) = 1;
        elseif ~isempty(labels(labels ~= 0)) && all(labels(labels ~= 0) == -1)
            gtConsolidated(i).labels(j) = -1;
        else
            gtConsolidated(i).labels(j) = 0;
        end
        prob = allProb(theseIdx);
        gtConsolidated(i).prob(j) = max(prob);
    end
    keepIdx = gtConsolidated(i).labels ~= 0;
    gtConsolidated(i).edges = gtConsolidated(i).edges(keepIdx,:);
    gtConsolidated(i).labels = gtConsolidated(i).labels(keepIdx);
    gtConsolidated(i).prob = gtConsolidated(i).prob(keepIdx);
end

end
