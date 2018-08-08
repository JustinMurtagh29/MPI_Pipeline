function [stats, debug] = gtReconstructionStats(skels, segIds, agglos, ...
    ov, opts)
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
%       opts: (Optional) struct
%           Struture with additional options.
%           'nodeDist': Double
%               If provided, then the reconstruction recall is also
%               evaluated using flight paths. A ground truth node is
%               considered "recalled", if there is any superagglo nodes
%               within the specified nodeDist. To run this the agglos have
%               to be in the superagglo format.
%               (Default: only "volume recall")
%           'voxelSize' [1x3] double
%               The voxel size in units of nodeDist.
%               (Default: skels{1}.scale)
%
% OUTPUT stats: table
%           Table with reconstruction statistics. All path lenghts are with
%           respect to skels{1}.scale.
%        debug: struct
%           Struct with some debugging outputs.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

defOpts = defaultOptions();
defOpts.voxelSize = skels{1}.scale;

if ~exist('opts', 'var') || isempty(opts)
    opts = defOpes;
else
    opts = Util.setUserOptions(defOpts, opts);
end

if isstruct(agglos)
    sagglos = agglos;
    agglos = Superagglos.getSegIds(agglos);
    isSuperagglo = true;
else
    isSuperagglo = false;
end

% get all nodes that overlap with
ovIds = cellfun(@(x)cell2mat(agglos(x(:,1))), ov, 'uni', 0);
recalledNodes_vol = cellfun(@(x, y) any(ismember(x, y), 2), segIds, ...
    ovIds, 'uni', 0);

if ~isempty(opts.nodeDist) && isSuperagglo
    getNodes = @(x)cell2mat(Superagglos.getNodes(x));
    recalledNodes_nodeDist = cellfun( ...
        @(skel, y)minDist(skel.getNodes(1), getNodes(sagglos(y(:, 1))), ...
        opts.voxelSize) < opts.nodeDist, skels, ov, 'uni', 0);
    hasNodeRec = true;
else
    recalledNodes_nodeDist = [];
    hasNodeRec = false;
end

ignoredNodes = cellfun(@(x)all(x < 0, 2), segIds, 'uni', 0);

recalledEdges_vol = cellfun(@(s, isRec) any(isRec(s.edges{1}), 2), ...
    skels, recalledNodes_vol, 'uni', 0);
ignEdges = cellfun(@(s, ign)any(ign(s.edges{1}), 2), skels, ...
    ignoredNodes, 'uni', 0);

% estimate mergers by distance to GT skeletons
isMerger = cell(length(ov), 1);
mergerNodes = cell(length(ov), 1);
for i = 1:length(ov)
    potMergeSeg = cellfun(@(x) ~ismember(x, segIds{i}), ...
        agglos(ov{i}(:, 1)), 'uni', 0);
    if isSuperagglo
        nodes = arrayfun(@Superagglos.getNodes, sagglos(ov{i}(:, 1)), ...
            'uni', 0);
    elseif ~isempty(opts.seg_com)
        nodes = cellfun(@(x)opts.seg_com(x, :), agglos(ov{i}(:, 1)), ...
            'uni', 0);
    else
        num_merger = nan(length(ov), 1);
        break;
    end
    
    d = cellfun(@(x)minDist(x, skels{i}.getNodes(1), opts.voxelSize), ...
        nodes, 'uni', 0);
    for j = 1:length(d)
        d{j}(~potMergeSeg{j}) = 0;
    end
    isMerger{i} = cellfun(@(x)any(x > opts.merge_dist), d);
    mergerNodes{i} = cellfun(@(x)x > opts.merge_dist, d, 'uni', 0);
    num_merger = cellfun(@(x)sum(x), isMerger);
end

% length of each edge
l = cellfun(@(x)edgeLength(x.nodes{1}(:, 1:3), x.edges{1}, opts.voxelSize), ...
    skels, 'uni', 0);

stats.pathLength = cellfun(@(x, ign)sum(x(~ign)), l, ignEdges);
stats.recalledPathLength_vol = cellfun(@(x, y, ign)sum(x(y & ~ign)), l, ...
    recalledEdges_vol, ignEdges);
stats.recall_volume = stats.recalledPathLength_vol ./ stats.pathLength;
if hasNodeRec
    recalledEdges_nodeDist = cellfun(@(s, isRec) any(isRec(s.edges{1}), 2), ...
        skels, recalledNodes_nodeDist, 'uni', 0);
    stats.recalledPathLength_nodeDist = cellfun( ...
        @(x, y, ign)sum(x(y & ~ign)), l, recalledEdges_nodeDist, ignEdges);
    stats.recall_nodeDist = stats.recalledPathLength_nodeDist ./ ...
        stats.pathLength;
end
stats.splits = cellfun(@(x)size(x, 1), ov);
stats.merger = num_merger;

stats = struct2table(stats);

if nargout > 1
    debug.recalledNodes_vol = recalledNodes_vol;
    debug.recalledNodes_nodeDist = recalledNodes_nodeDist;
    debug.ignoredNodes = ignoredNodes;
    debug.mergerNodes = mergerNodes;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = minDist(x, y, voxelSize)
% INPUT x: [Nx3] double
%           The ground truth skeleton nodes
%       y: [Nx3] double
%           The nodes of an agglo overlapping with the skeleton.

if ~all(voxelSize == 1)
    x = double(x) .* voxelSize(:)';
    y = double(y) .* voxelSize(:)';
end
d = min(pdist2(x, y), [], 2);
end

function l = edgeLength(nodes, edges, vxSize)
% the length for the single edges

l = nodes(edges(:, 1), :) - nodes(edges(:, 2), :);
l = l .* vxSize(:)';
l = sqrt(sum(l.^2, 2));

end

function opts = defaultOptions()
opts.nodeDist = [];
opts.voxelSize = [];
opts.seg_com = [];
opts.merge_dist = 1000;
end