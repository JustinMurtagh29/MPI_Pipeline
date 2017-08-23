function [agglo, seed2agglo] = distRestrictedAgglo( seeds, edges, ...
    contProb, probT, dist, combineOverlapping )
%DISTRESTRICTEDAGGLO Agglomeration via contuinty probs with a distance
%restriction.
% INPUT seeds: [Nx1] int
%           The segment ids of the seeding segments.
%       edgesOrNeighbors: [Nx2] int or [Nx1] cell
%           The graph edge list or directly the neighbors list.
%       contProb: [Nx1] double
%           The continuity probabilites (= edge weights) for the
%           corresponding edge.
%           (This is only required if the edges are supplied as second
%           input).
%       probT: double
%           Lower threshold on the contProb to consider segments as
%           connected.
%           (This is only required if the edges are supplied as second
%           input).
%       dist: int
%           The maximal distance on the graph to the seed segments.
%       combineOverlapping: (Optional) logical
%           Flag to combine overlapping agglos. The second output in this
%           case contains the indices from seed to agglo.
%           (Default: false)
% OUTPUT agglo: [Nx1] cell
%           Cell array containing the output agglos. Depending on the
%           overlapMode this has the same length as seeds ('None') or is
%           shorter ('combine').
%        seed2agglo: [Nx1] int
%           Linear index of the output agglo for the corresponding input
%           seed (only for mode 'getIdx').
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if iscell(edges)
    nIds = edges;
else
    % get edges above threshold
    edgesT = edges(contProb > probT, :);

    % get the neighbor list
    nIds = Graph.edges2Neighbors(edgesT);
end

% do the agglomeration based on the neighbors
agglo = num2cell(seeds);
for i = 1:dist
    agglo = cellfun(@(x)unique([x; cell2mat(nIds(x))], 'stable'), ...
        agglo, 'uni', 0);
end

if exist('combineOverlapping', 'var') && combineOverlapping
    [agglo, seed2agglo] = Seg.Global.combineEClasses(agglo, false);
else
    seed2agglo = [];
end

end

