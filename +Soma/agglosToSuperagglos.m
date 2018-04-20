function somata = agglosToSuperagglos( p, outFile )
%AGGLOSTOSUPERAGGLOS Convert to some agglos to the superagglo format.
% INPUT p: struct
%           Segmentation parameter struct.
%           Soma agglos are loaded from p.agglo.somaFile
%       outFile: (Optional) string
%           Path to output file where the superagglos are stored.
% OUTPUT somaSagglos: [Nx1] struct
%           Soma agglos in the superagglo format.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

info = Util.runInfo(false);

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
    id2Idx = zeros(maxId, 1);
    id2Idx(thisIds) = 1:length(thisIds);
    thisEdgesAgglo = id2Idx(thisEdges);
    thisEdgesAgglo = unique(sort(thisEdgesAgglo, 2), 'rows');
    somata(i, 1).nodes = [pt(thisIds,:), thisIds];
    somata(i, 1).edges = double(thisEdgesAgglo);
    cc = Graph.findConnectedComponents(somata(i, 1).edges);
    if length(cc) > 1
        % add random edges between components
        conn = [repelem(cc{1}, length(cc) - 1)', ...
            cellfun(@(x)x(1), cc(2:end))];
        somata(i, 1).edges = cat(1, somata(i, 1).edges, conn);
    end
end

if exist('outFile', 'var') && ischar(outFile)
    Util.log('Saving output to %s.', outFile);
    Util.save(outFile, somata, centerSomaIdx, info);
end

end

