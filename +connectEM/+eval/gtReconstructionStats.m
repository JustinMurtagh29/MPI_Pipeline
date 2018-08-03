function [stats, debug] = gtReconstructionStats(skels, segIds, agglos, ov)
%GTRECONSTRUCTIONSTATS Statistics about the reconstruction of the ground
% truth axons.
%
% INPUT skels: [Nx1] cell
%           The ground truth skeletons.
%       segIds: [Nx1] cell
%           Segment ids for the nodes of the ground truth skeletons.
%           (see also second output of connectEM.eval.getNewAxonGT)
%       agglos: [Nx1] cell or struct
%           The reconstructed agglomerates in the agglo or superagglo
%           format.
%       ov: [Nx1] cell
%           Cell with linear indices of agglos that overlap with the
%           corresponding gtSegIds.
%           (see output of connectEM.eval.getNewAxonGTAggloOverlap)
%
% OUTPUT stats: struct
%           Struct with reconstruction statistics.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if isstruct(agglos)
    agglos = Superagglos.getSegIds(agglos);
end

% get all nodes that overlap with
ovIds = cellfun(@(x)cell2mat(agglos(x(:,1))), ov, 'uni', 0);
recalledNodes = cellfun(@(x, y) ismember(x, y), segIds, ovIds, 'uni', 0);
ignoredNodes = cellfun(@(x)x < 0, segIds, 'uni', 0);

recalledEdges = cellfun(@(s, isRec) any(isRec(s.edges{1}), 2), ...
    skels, recalledNodes, 'uni', 0);
ignEdges = cellfun(@(s, ign)any(ign(s.edges{1}), 2), skels, ...
    ignoredNodes, 'uni', 0);

% length of each edge
l = cellfun(@(x)edgeLength(x.nodes{1}(:, 1:3), x.edges{1}, x.scale), ...
    skels, 'uni', 0);

stats.pathLength = cellfun(@(x, ign)sum(x(~ign)), l, ignEdges);
stats.recalledPathLength = cellfun(@(x, y, ign)sum(x(y & ~ign)), l, ...
    recalledEdges, ignEdges);
stats.recall = stats.recalledPathLength ./ stats.pathLength;

debug.recalledNodes = recalledNodes;
debug.ignoredNodes = ignoredNodes;

end

function l = edgeLength(nodes, edges, vxSize)
% the length for the single edges

l = nodes(edges(:, 1), :) - nodes(edges(:, 2), :);
l = l .* vxSize(:)';
l = sqrt(sum(l.^2, 2));

end