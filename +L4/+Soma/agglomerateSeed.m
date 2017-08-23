function [ids, cEdges] = agglomerateSeed( seedId, edges, mergeP, t )
%AGGLOMERATESEED Agglomerate from a seed segment.
% INPUT seedId: int
%           Segmentation id of the starting segment.
%       edges: [Nx2] int
%           Edge list.
%       mergeP: [Nx1] double
%           Merge probability for the corresponding edges.
%       t: double
%           Lower probability threshold for merging.
% OUTPUT ids: [Nx1] int
%           Ids that are connected to the seedId.
%        eEdges: [Nx2] int
%           The edges that connect the ids.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

cEdges = edges(mergeP > t, :);
cc = Graph.findConnectedComponents(cEdges, false, false);
ids = cc{cellfun(@(x)any(x == seedId), cc)};

end

