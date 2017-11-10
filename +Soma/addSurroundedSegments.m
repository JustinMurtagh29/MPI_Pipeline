function segIds = addSurroundedSegments( agglos, edges, borderSize, ...
    fracT, corrEdges, iter )
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
%       corr: (Optional) [Nx2] int
%           Edge list of correspondences. If supplied then the borderSize
%           of segments is considered w.r.t. all borders of corresponding
%           segments.
%           (Default: no correspondences)
%       iter: (Optional) int
%           Number of iterations that the inclusion of surrounding segments
%           is performed. If there are no segments added in one iterations
%           the the iterations are automatically stopped.
%           (Note that this can have an effect, e.g. when
%           there is one large and one small missed segment in the soma 
%           that are neighboring and the small segment has a large surface
%           area with the large one).
%           (Default: 1)
% OUTPUT segIds: [Nx1] cell
%           Cell array of segment ids that are surrounded by the
%           corresponding agglo.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if exist('corrEdges', 'var') && ~isempty(corrEdges)
    % concatenate correspondences in neighbor lists
    eLUT = Seg.Global.EdgeLookup(edges, true);
    eClasses = Graph.findConnectedComponents(corrEdges, false, false);
    for i = 1:length(eClasses)
        edges(eLUT.edgeIdx(eClasses{i}, false)) = eClasses{i}(1);
    end
end

if ~exist('iter', 'var') || isempty(iter)
    iter = 1;
end

[nIds, nIdx] = Graph.edges2Neighbors(edges);

segIds = cell(length(agglos), 1);
for i = 1:length(agglos)
    for it = 1:iter
        isInAgglo = false(max(edges(:)), 1);
        isInAgglo(agglos{i}) = true;
        bordOnAgglo = cellfun( ...
            @(x, y)sum(borderSize(x(isInAgglo(y))))/sum(borderSize(x)), ...
            nIdx, nIds);
        toAddIds = setdiff(find(bordOnAgglo >= fracT), agglos{i});
        segIds{i} = cat(1, segIds{i}, toAddIds);
        agglos{i} = cat(1, agglos{i}(:), segIds{i}(:));
        if isempty(toAddIds)
            break;
        end
    end
end

end

