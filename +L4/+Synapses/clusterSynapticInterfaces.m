function synapses = clusterSynapticInterfaces( edges, synScores, ...
    synIdx, contProb, probT, dist )
%CLUSTERSYNAPTICINTERFACES Synapse clustering by local agglomeration of
% pre- and postsynaptic processes.
% INPUT edges: [Nx2] int
%           Global edge list.
%       synScores: [Nx2] float
%           Synapse scores for the corresponding edges.
%       synIdx: [Nx1] logical
%           The detected synaptic edges.
%       contProb: [Nx1] float
%           Merge continuity probability for the correspoding edge.
%       probT: float
%           Threhold for local bouton agglomeration.
%       dist: int
%           Maximal distance in number of edges that is agglomeration
%           starting from a presynaptic segment.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

% get pre- and postsyn ids
Util.log('Getting synaptic interfaces.');
synEdges = edges(synIdx, :);
synScores = synScores(synIdx,:);
idx1 = synScores(:,1) >= synScores(:,2);
presynId = synEdges([idx1, ~idx1]);
postsynId = synEdges([find(idx1) + size(synEdges, 1); find(~idx1)]);
edgeIdx = find(synIdx);
edgeIdx = edgeIdx([find(idx1); find(~idx1)]);

% edges and neighbors above probT
Util.log('Preparing agglomeration graph.');
edgesT = edges(contProb > probT, :);
nIds = Graph.edges2Neighbors(edgesT);

% do the agglomeration
Util.log('Running presynaptic agglomeration.');
[~, preAggloIdx] = L4.Agglo.distRestrictedAgglo(presynId, nIds, [], ...
    [], dist, true);
Util.log('Running postsynaptic agglomeration.');
[~, postAggloIdx] = L4.Agglo.distRestrictedAgglo(postsynId, nIds, [], ...
    [], dist, true);

% combine synaptic interfaces to synapses based on the overlapping local
% agglomerations

Util.log('Combining interfaces between same pre- and postsynaptic agglo.');
[~,~,synIdx] = unique([preAggloIdx, postAggloIdx], 'rows');
synapses = accumarray(synIdx, edgeIdx, [], @(x){x});

end

