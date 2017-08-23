function synIdx = getSynapsesForAgglos(sAgglos, synIdx, edges)
%GETSYNAPSESFORAGGLOS Get the edge indices for synapses along agglos.
% INPUT sAgglos: [Nx1] struct array
%           The superagglos struct array.
%       synIds: [Nx1] logical
%           Logical indices for synaptic edges w.r.t. global edge list.
%       edges: [Nx2] int
%           The global edge list.
% OUTPUT synIdx: [Nx1] cell
%           The linear indices w.r.t. the global edges that are synaptic
%           for the corresponding superagglo.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

edgeLookup = Seg.Global.getEdgeLookupTable(edges);
synEdgeLookup = cellfun(@(x)x(synIdx(x)), edgeLookup, 'uni', 0);

synIdx = cell(length(sAgglos), 1);
for i = 1:length(sAgglos)
    synIdx{i} = cell2mat(synEdgeLookup(getSegIds(sAgglos(i))));
end
end

function ids = getSegIds(sAgglo)
    ids = sAgglo.nodes((sAgglo.nodes(:,4) > 0), 4);
end
