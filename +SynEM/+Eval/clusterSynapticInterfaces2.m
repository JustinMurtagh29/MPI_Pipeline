function T = clusterSynapticInterfaces2( coms, cutoff, scale, chunkSize )
%CLUSTERSYNAPTICINTERFACES Cluster interface coms based on single linkage.
% INPUT coms: [Nx3] float
%           Global coms of the synaptic interfaces.
%       cutoff: double
%           Distance cutoff.
%       scale: (Optional) [1x3] double
%           Scaling factor for each dimension. cutoff should be specified
%           in units of scale.
% OUTPUT T: [Nx1] int
%           Cluster index for each input com.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if exist('scale', 'var') && ~isempty(scale)
    coms = bsxfun(@times, single(coms), scale(:)');
else
    coms = single(coms);
end

if ~exist('chunkSize', 'var') || isempty(chunkSize)
    chunkSize = 1e5;
end

l = size(coms, 1);
adj = false(l, l);
numChunks = ceil(l / chunkSize);
for i = 1:numChunks
    for j = 1:numChunks
        from_i = (i-1)*chunkSize + 1;
        to_i = min(chunkSize*i, l);
        from_j = (j-1)*chunkSize + 1;
        to_j = min(chunkSize*j, l);
        adj(from_i:to_i, from_j:to_j) = ...
            pdist2(coms(from_i:to_i,:), coms(from_j:to_j, :)) < cutoff;
    end
end

adj = sparse(adj);
adj = adj + adj';
[~, T] = Graph.findConnectedComponents(adj, false);
T = T(:);

end

