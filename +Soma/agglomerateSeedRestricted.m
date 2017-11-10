function [ids, cEdges] = agglomerateSeedRestricted( seedIds, seedC, ...
    edges, segCom, borderSize, mergeP, tp, ts, borderCom, r, ...
    voxelSize )
%AGGLOMERATESEEDRESTRICTED Agglomerate around a seed id restricted to a
%given distance.
% INPUT seedIds: [Nx1]
%           Segmentation ids of the starting segments.
%       seedC: [1x3] int/float
%           Global coordinates of the seed id.
%       edges: [Nx2] int
%           Edge list.
%       segCom: [3xN] double
%           Global segment coms.
%       borderSize: [Nx1] double
%           The size for each border.
%       mergeP: [Nx1] double
%           Merge probability for the corresponding edges.
%       tp: double
%           Lower probability threshold for agglomeration.
%       ts: double
%           Lower border size threshold for agglomeration.
%       borderCoM: [Nx3] int/float
%           Coms of the corresponding edges. Should have the same units as
%           seedC.
%       r: double
%           Maximal distance of borders to seedC (in units of
%           seedC/borderCom).
%       voxelSize: (Optional) [1x3] double
%           Scaling factor of the global coordinates.
%           (Default: [1, 1, 1])
% OUTPUT ids: [Nx1] int
%           Ids that are connected to the seedId.
%        cEdges: [Nx2] int
%           The edges that connect the ids.
% Author: Benedikt Staffler, modified by Robin Hesse

%% preprocess coordinates
seedC = single(seedC(:)');
if exist('voxelSize', 'var') && ~isempty(voxelSize)
    voxelSize = voxelSize(:)';
    seedC = seedC.*voxelSize;
    borderComS = uint64(borderCom);
    borderComS = bsxfun(@times, single(borderComS), voxelSize);
end

idx = find(pdist2(borderComS, seedC) < r);

edges = edges(idx, :);
borderSize = borderSize(idx);
mergeP = mergeP(idx);

%% remove distant cube correspondences

edgePoints1 = segCom(:,edges(:,1))';
edgePoints2 = segCom(:,edges(:,2))';
edgePoints1 = bsxfun(@times, single(edgePoints1), voxelSize);
edgePoints2 = bsxfun(@times, single(edgePoints2), voxelSize);

idx = find(pdist2(edgePoints1, seedC) < r | pdist2(edgePoints2, seedC) < r );

[ids, cEdges] = Soma.agglomerateSeed(seedIds, edges(idx, :),borderSize(idx),...
	mergeP(idx), tp, ts);


end

