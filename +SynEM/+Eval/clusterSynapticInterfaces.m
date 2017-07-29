function T = clusterSynapticInterfaces( coms, cutoff, scale )
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

numPoints = size(coms, 1);
if exist('scale', 'var') && ~isempty(scale)
    coms = bsxfun(@times, single(coms), scale(:)');
else
    coms = single(coms);
end

try
    T = clusterdata(coms, 'linkage', 'single', ...
            'criterion', 'distance', 'cutoff', cutoff);
catch
    % handle case of too much memory requirement for direct clusterdata
    adj = cell(size(coms, 1), 1);
    for i = 1:(numPoints - 1)
        adj{i} = i + find(pdist2(coms(i,:), coms(i+1:end,:)) < cutoff)';
    end
    l = cellfun(@length, adj);
    adj = cell2mat(adj);
    adjM = sparse(repelem(1:numPoints, l), adj, 1, numPoints, numPoints);
    [~, T] = Graph.findConnectedComponents(adjM + adjM', false);
    T = T(:);
end


end

