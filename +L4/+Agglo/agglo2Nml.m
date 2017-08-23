function [ skel ] = agglo2Nml( agglos, segmentCom )
%AGGLO2NML Convert agglomerations to a nml using MST of centroids.
% INPUT agglos: [Nx1] int or [Nx1] cell
%           Segment ids of the agglomeration segments or cell array of
%           segment ids. Each cell will be saved to a separate tree.
%       segmentPoints: [Nx3] int
%           Segment center of mass/points that will be used as the node
%           location for an agglomeration id.
% OUTPUT skel: skeleton object
%           Skeleton object containing each agglomeration as a single tree.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~iscell(agglos)
    agglos = {agglos};
end

skel = Skeleton.fromMST( ...
    cellfun(@(ids) {segmentCom(ids, :)}, agglos), [11.24, 11.24, 28]);
skel = Skeleton.setParams4Dataset(skel, 'ex145_ROI2017');

end

