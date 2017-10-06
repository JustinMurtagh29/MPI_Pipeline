function [ skels, segIds ] = getNewAxonGT(p)
%GETNEWAXONGT Load the axons in the evaluationData/new_axon_gt/ folder.
% INPUT p: (Optional) struct
%           Segmentation parameter struct for ROI2017.
%           (Default: Gaba.getSegParameters('ex145_ROI2017'))
% OUTPUT skels: [10x1] cell
%           Cell array with skeleton objects for each new_axon_gt nml file.
%        segIds: [Nx1] cell
%           Segment ids for the skeletons in skels. If only calculated if
%           nargout > 1.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

thisPath = mfilename('fullpath');
ax_gt_path = fileparts(fileparts(thisPath));
ax_gt_path = fullfile(ax_gt_path, 'evaluationData', 'new_axon_gt');

if ~exist('p', 'var') || isempty(p)
    p = Gaba.getSegParameters('ex145_ROI2017');
end

skels = skeleton.loadSkelCollection(ax_gt_path, [], true);
if nargout > 1
    segIds = cellfun(@(x)Skeleton.getSegmentIdsOfSkel(p, x), skels, ...
        'uni', 0);
    % each nml has exactly one tree
    segIds = cellfun(@(x)x{1}, segIds, 'uni', 0);
end
end

