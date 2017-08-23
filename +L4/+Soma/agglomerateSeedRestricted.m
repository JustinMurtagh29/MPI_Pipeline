function [ids, cEdges] = agglomerateSeedRestricted( seedId, seedC, ...
    edges, mergeP, t, borderCom, r, voxelSize )
%AGGLOMERATESEEDRESTRICTED Agglomerate around a seed id restricted to a
%given distance.
% INPUT seedId: int
%           Segmentation id of the starting segment.
%       seedC: [1x3] int/float
%           Global coordinates of the seed id.
%       edges: [Nx2] int
%           Edge list.
%       mergeP: [Nx1] double
%           Merge probability for the corresponding edges.
%       t: double
%           Lower probability threshold for merging.
%       borderCoM: [Nx3] int/float
%           Coms of the corresponding edges. Should have the same units as
%           seedC.
%       r: double
%           Maximal distance of borders to seedC (in units of
%           seedC/borderCom).
%       voxelSize: (Optional) [1x3] double
%           Scaling factor of the global coordinates.
%           (Default: [1, 1, 1])
%        eEdges: [Nx2] int
%           The edges that connect the ids.
% OUTPUT ids: [Nx1] int
%           Ids that are connected to the seedId.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

% preprocess coordinates
seedC = single(seedC(:)');
if exist('voxelSize', 'var') && ~isempty(voxelSize)
    voxelSize = voxelSize(:)';
    seedC = seedC.*voxelSize;
    borderCom = bsxfun(@times, single(borderCom), voxelSize);
end

idx = pdist2(borderCom, seedC) < r;
[ids, cEdges] = L4.Soma.agglomerateSeed(seedId, edges(idx, :), ...
    mergeP(idx), t);


end

