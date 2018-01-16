function sagglo = mergeMinDist( sagglo1, sagglo2, scale )
%MERGEMINDIST Merge two superagglos at the nodes with minimal distance.
% INPUT sagglo1: struct
%           First superagglo.
%       sagglo2: struct
%           Second superagglo.
%       sacles: (Optional) [1x3] double
%           Node scale used for distance calculation.
%           (Default: [11.24, 11.24, 28])
% OUTPUT sagglo: struct
%           The merged superagglo.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if nargin < 3
    scale = [11.24, 11.24, 28];
end
scale = scale(:)';

% cat fields from old agglos
nodes = cat(1, sagglo1.nodes, sagglo2.nodes);
offset = size(sagglo1.nodes, 1);
edges = cat(1, sagglo1.edges, sagglo2.edges + offset);
if isfield(sagglo1, 'endings') && isfield(sagglo2, 'endings')
    endings = cat(1, sagglo1.endings, sagglo2.endings + offset);
end
if isfield(sagglo1, 'solvedChiasma') && isfield(sagglo2, 'solvedChiasma')
    solvedChiasma = cat(1, sagglo1.solvedChiasma, sagglo2.solvedChiasma);
end

% introduce a new node with shortest distance
d = pdist2(bsxfun(@times, sagglo1.nodes, scale), ...
           bsxfun(@times, sagglo1.nodes, scale));
[~, idx] = min(d);
[x, y] = ind2sub(size(d), idx);

% add new edge between shortest nodes
edges(end+1, :) = [x, y + offset];

% create new agglo
sagglo.nodes = nodes;
sagglo.edges = edges;
sagglo.endings = endings;
sagglo.solvedChiasma = solvedChiasma;

end

