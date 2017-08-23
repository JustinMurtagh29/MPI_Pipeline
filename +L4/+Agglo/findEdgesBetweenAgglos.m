function [c, edgeIdx] = findEdgesBetweenAgglos( agglos, edges )
%FINDEDGESBETWEENAGGLOS
% INPUT agglos: [Nx1] cell
%           Each cell contains the ids of one agglomeration.
%       edges: [Nx2] int
%           Segment adjacency graph edges.
% OUTPUT c: [Nx2] int
%           The agglos that have edges between them.
%        edgeIdx: [Nx1] cell
%           The linear edge indices of edges for the corresponding agglos.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

% edge list wrt agglos
lut = L4.Agglo.buildLUT(max(edges(:)), agglos);
edgesA = lut(edges);

% edge indices between different agglos
idx = all(edgesA, 2) & (edgesA(:,1) ~= edgesA(:,2));

edgesA = edgesA(idx, :);
edgesA = sort(edgesA, 2);
[c, ~, ic] = unique(edgesA, 'rows');
edgeIdx = accumarray(ic, find(idx), [], @(x){x});

end

