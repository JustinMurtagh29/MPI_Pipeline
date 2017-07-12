function [graph, segmentMeta, borderMeta, globalSegmentPCA] = loadAllSegmentationData(p);

    graph = load([p.saveFolder 'graphNew.mat'], 'edges', 'prob', 'borderIdx');
    [graph.neighbours, neighboursIdx] = Graph.edges2Neighbors(graph.edges);
    graph.neighProb = cellfun(@(x)graph.prob(x), neighboursIdx, 'uni', 0);
    graph.neighBorderIdx = cellfun(@(x)graph.borderIdx(x), neighboursIdx, 'uni', 0);

    if nargout > 1
        segmentMeta = load([p.saveFolder 'segmentMeta.mat.20170523'], '-mat');
        segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);
    end

    if nargout > 2
        borderMeta = load([p.saveFolder 'globalBorder.mat'], 'borderSize', 'borderCoM');
    end

    if nargout > 3
        globalSegmentPCA = load([p.saveFolder 'globalSegmentPCA.mat'], 'covMat');
    end

end

