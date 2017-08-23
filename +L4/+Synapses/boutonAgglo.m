function [bAgglo, postsynId, edgeIdx] = boutonAgglo(edges, synScores, synIdx, ...
    mergeProb, mergeT, dist, combineOverlapping)
%BOUTONAGGLO Location agglomeration in presynptic segments.
% INPUT edges: [Nx2] int
%           Global edge list.
%       synScores: [Nx2] float
%           Synapse scores for the corresponding edges.
%       synIdx: [Nx1] logical
%           The detected synaptic edges.
%       mergeProb: [Nx1] float
%           Merge continuity probability for the correspoding edge.
%       mergeT: float
%           Threhold for local bouton agglomeration.
%       dist: int
%           Maximal distance in number of edges that is agglomeration
%           starting from a presynaptic segment.
%       combineOverlapping: (Optional) logical
%           Flag indicating whether overlapping boutons should be combined.
%           Boutons are considered overlaping if they have a single segment
%           in common.
% OUTPUT bAgglo: [Nx1] cell
%           Cell array of bouton agglos. The first id in bAgglo is the
%           detected presynaptic segment.
%        postsynId: [Nx1] int
%           The index of the corresponding postsynapic segment for a
%           bouton. In case combineOverlapping is true then postsynId is a
%           cell array with the postsynaptic segment ids for each bouton.
%        edgeIdx: [Nx1] int
%           The global linear index of the edge connecting the
%           corresponding bouton with the corresponding postsynId.
%           In case combineOverlapping is true then edgeIdx is a
%           cell array with the synaptic edge indices for each bouton.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

% get presyn ids
synEdges = edges(synIdx, :);
synScores = synScores(synIdx,:);
idx1 = synScores(:,1) >= synScores(:,2);
presynId = synEdges([idx1, ~idx1]);
postsynId = synEdges([find(idx1) + size(synEdges, 1); find(~idx1)]);
edgeIdx = find(synIdx);
edgeIdx = edgeIdx([find(idx1); find(~idx1)]);

% do a local agglomeration
edgesT = edges(mergeProb > mergeT, :);

% get neighbors
nIds = Graph.edges2Neighbors(edgesT);

% do the agglomeration based on the neighbors
bAgglo = num2cell(presynId);
for i = 1:dist
    bAgglo = cellfun(@(x)unique([x; cell2mat(nIds(x))], 'stable'), ...
        bAgglo, 'uni', 0);
end

if exist('combineOverlapping', 'var') && combineOverlapping
    [bAgglo, idx] = Seg.Global.combineEClasses(bAgglo, false);
    postsynId = accumarray(idx, postsynId, [], @(x){x});
    edgeIdx = accumarray(idx, edgeIdx, [], @(x){x});
end

end
