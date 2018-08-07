function nodes = getNodes(sagglos)
%GETNODES Return the nodes of the input superagglos.
%
% INPUT sagglos: [Nx1] struct
%           The superagglos
%
% OUTPUT nodes: [Nx1] cell
%           The nodes as [Nx3] int for each input superagglo.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

nodes = arrayfun(@(x)x.nodes(:, 1:3), sagglos, 'uni', 0);

end

