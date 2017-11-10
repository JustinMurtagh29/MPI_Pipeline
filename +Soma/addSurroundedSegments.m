function segIds = addSurroundedSegments( agglos, edges, borderSize, fracT )
%ADDSURROUNDEDSEGMENTS Find segments for which a fractions of the
% borderArea (currently in number of voxels) is onto an agglo.
% INPUT agglos: [Nx1] cell
%           Cell array of segment ids.
%       edges: [Nx2] int
%           Global edge list.
%       borderSize: [Nx1] double
%           The border size for each edge.
%       fracT: double
%           Fraction of the total border of a segment that is covered by
%           the agglo to consider the agglo surrounded.
% OUTPUT segIds: [Nx1] cell
%           Cell array of segment ids that are surrounded by the
%           corresponding agglo.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

neighbors = Graph.edges2Neighbors(edges);
segIds = cell(length(agglos), 1);
for i = 1:length(agglos)
    isInAgglo = false(max(edges(:)), 1);
    isInAgglo(agglos{i}) = true;
    bordOnAgglo = cellfun( ...
        @(x)sum(borderSize(x(isInAgglo(x))))/sum(borderSize(x)), neighbors);
    segIds{i} = setdiff(find(bordOnAgglo >= fracT), agglos{i});
end

end

