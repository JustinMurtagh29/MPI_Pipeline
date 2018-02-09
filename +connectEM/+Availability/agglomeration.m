function newAgglo = agglomeration (existingAgglo, borderMeta, graph, segmentMeta)
    % set all existing Agglo pieces to true
    segmentMeta.axonProb(existingAgglo) = 1;
    segmentMeta.voxelSize(existingAgglo) = Inf;
    
    % greenlight the NaNs (= Correspondences)
    borderMeta.borderSize(end+1) = Inf
    borderMeta.borderCoM(end+1, :) = [NaN, NaN, NaN];
    graph.borderIdx(graph.borderIdx == NaN) = length(borderMeta.borderSize);
    
    
    %classic agglomeration
    edges = graph.prob > 0.7;
    edges = edges & borderMeta.borderSize(graph.borderIdx) > 30;
    edges = edges & all(segmentMeta.axonProb(graph.edges) > 0.3, 2);
    edges = edges & all(segmentMeta.voxelCount(graph.edges) > 10, 2);
    
    %distance threshold
    scalize = @(v)bsxfun(@times, v, [11.24, 11.24, 26]);
    hotEdges = borderMeta.borderCoM(graph.borderIdx(edges));
    filterNan = @(x)x(~isnan(x))
    surface = borderMeta.borderCoM(filterNan(graph.borderIdx(any(ismember(graph.edges,existingAgglo,'rows'),2))));
    edges(edges) = any(min(pdist2(scalize(hotEdges), scalize(surface)))<2000|isnan(hotEdges))  ;
    
    % do the agglomeration
    C = Graph.findConnectedComponents([graph.edges(edges); existingAgglo(2:end), existingAgglo(1:end-1)]);
    findIdx = find(cellmat(@(x)ismember(x,existingAgglo(1))));
    newAgglo = C{findIdx};