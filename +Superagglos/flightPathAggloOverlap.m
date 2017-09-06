function ov = flightPathAggloOverlap( p, fp, agglos, nodeEv )
%FLIGHTPATHAGGLOOVERLAP Get the overlap of flight paths with agglos.
% INPUT p: struct
%           Segmentation parameter struct.
%       fp: [Nx1] struct
%           Flight paths in the super-agglo format (only nodes are actually
%           needed).
%           (see also Superagglos.getFlightPath)
%       agglos: [Nx1] cell or [Nx1] struct
%           The agglos/superagglos for which the overlap with the flight
%           paths are calculated.
%       nodeEv: (Optional) double
%           Minimal required node evidence for an agglo to be picked up.
% OUTPUT ov: [Nx1] cell
%           Cell array containing the linear indices of the agglos that are
%           overlapping with the corresponding flight path.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if isstruct(agglos)
    agglos = Superagglos.getSegIds(agglos);
end

aggloLUT = L4.Agglo.buildLUT(agglos);
fpIds = getFPNodeSegIds(p, fp);

% get nodes in agglos
fpAggloIdx = cellfun(@(x)aggloLUT(x(x > 0)), fpIds, 'uni', 0);
fpAggloIdx = cellfun(@(x)x(x > 0), fpAggloIdx, 'uni', 0);

% node evidence
[c, nc] = cellfun(@Util.ucounts, fpAggloIdx, 'uni', 0);

% get overlapping agglos
ov = cellfun(@(x, y) x(y >= nodeEv), c, nc, 'uni', 0);

end

function ids = getFPNodeSegIds(p, fp)

nodes = {fp.nodes}';
l = cellfun(@(x)size(x, 1), nodes);
nodes = cell2mat(nodes);
ids = zeros(size(nodes, 1), 27);
cubeIdx = Skeleton.findSegCubeIdxOfNodes(nodes, p);
[group, cubeIdx] = Util.group2Cell(cubeIdx);
for i = 1:length(cubeIdx)
    this_bbox = Util.addBorder(p.local(cubeIdx(i)).bboxSmall, [1, 1, 1]);
    seg = Seg.IO.loadSeg(p, this_bbox);
    curNodes = Util.indConversion(this_bbox, nodes(group(i),:));
    curNodes = bsxfun(@plus, curNodes, [0, Util.lneigh26(size(this_bbox))]);
    ids(group(i)) = seg(curNodes);
end
ids = mat2cell(ids, l, 27);

end