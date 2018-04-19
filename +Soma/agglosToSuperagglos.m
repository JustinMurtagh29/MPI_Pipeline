function somata = agglosToSuperagglos( p, outFile, probT )
%AGGLOSTOSUPERAGGLOS Convert to some agglos to the superagglo format.
% INPUT p: struct
%           Segmentation parameter struct.
%           Soma agglos are loaded from p.agglo.somaFile
%       outFile: (Optional) string
%           Path to output file where the superagglos are stored.
%       probT: (Optional) double
%           Continuity probability threshold used during agglomeration that
%           is used for the superagglo nodes as far as possible.
%           (Default: all edges between the agglo segments are used).
% OUTPUT somaSagglos: [Nx1] struct
%           Soma agglos in the superagglo format.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

info = Util.runInfo(false);

if ~exist('probT', 'var') || isempty(probT)
    probT = 0.98;
end

Util.log('Loading soma agglos from %s.', p.agglo.somaFile);
m = load(p.agglo.somaFile);
shAgglos = cellfun(@sort, m.somaAgglos(:,1), 'uni', 0);
centerSomaIdx = m.centerSomaIdx;

Util.log('Loading supervoxel graph.');
[graph, segmentMeta] = Seg.IO.loadGraph(p, false);
pt = segmentMeta.point;
maxId = max(graph.edges(:));

somata = struct();
Util.log('Converting soma agglos to superagglo format.');
for i = 1:length(shAgglos)
    thisIds = shAgglos{i};
    lut = L4.Agglo.buildLUT({thisIds}, maxId);
    idx = all(lut(graph.edges), 2);
    thisEdges = graph.edges(idx, :);
    aboveProb = graph.prob(idx) > probT;
    thisEdgesAboveProb = thisEdges(aboveProb, :);
    
    % for unconnected nodes use all svg edges
    toAddIds = setdiff(thisIds, thisEdgesAboveProb(:));
    toAdd = any(ismember(thisEdges, toAddIds), 2);
    thisAggloEdges = thisEdges(aboveProb | toAdd, :);
    
    id2Idx = zeros(maxId, 1);
    id2Idx(thisIds) = 1:length(thisIds);
    thisEdgesAgglo = id2Idx(thisAggloEdges);
    thisEdgesAgglo = unique(sort(thisEdgesAgglo, 2), 'rows');
    somata(i, 1).nodes = [pt(thisIds,:), thisIds];
    somata(i, 1).edges = double(thisEdgesAgglo);
end

if exist('outFile', 'var') && ischar(outFile)
    Util.log('Saving output to %s.', outFile);
    Util.save(outFile, somata, centerSomaIdx, info);
end

end

