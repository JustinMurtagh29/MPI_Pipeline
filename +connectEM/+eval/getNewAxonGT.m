function [ skels, segIds ] = getNewAxonGT(p, nhood)
%GETNEWAXONGT Load the axons in the evaluationData/new_axon_gt/ folder and
% transform them to the ROI2017 alignment.
%
% INPUT p: (Optional) struct
%           Segmentation parameter struct for ROI2017.
%           (Default: Gaba.getSegParameters('ex145_ROI2017'))
%       nhood: (Optional) int
%           Neighborhood for the extraction of segment ids.
%           (see also Skeleton.getSegmentIdsOfNodes)
%           (Default: 0)
% OUTPUT skels: [10x1] cell
%           Cell array with skeleton objects for each new_axon_gt nml file.
%        segIds: [Nx1] cell
%           Segment ids for the skeletons in skels. If only calculated if
%           nargout > 1.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('nhood', 'var') || isempty(nhood)
    nhood = 0;
end

thisPath = mfilename('fullpath');
ax_gt_path = fileparts(fileparts(thisPath));
ax_gt_path = fullfile(ax_gt_path, 'evaluationData', 'new_axon_gt');

if (~exist('p', 'var') || isempty(p)) && nargout > 1
    p = Gaba.getSegParameters('ex145_ROI2017');
end

skels = skeleton.loadSkelCollection(ax_gt_path, [], true);
skels = convertToROI2017Alignment(skels);

if nargout > 1
    segIds = cellfun( ...
        @(x)Skeleton.getSegmentIdsOfNodes(p, x.nodes{1}(:, 1:3), nhood), ...
        skels, 'uni', 0);
end
end

function skels = convertToROI2017Alignment(skels)
% code taken from connectEM.evaluateAggloMeta
for i = 1:10
    skels{i}.nodes{1}(:,1:3) = bsxfun(@minus, skels{i}.nodes{1}(:,1:3), ...
        [1195, 1515, 115] - (129 - [25 25 10]));
    skels{i}.nodes{1}(skels{i}.nodes{1} <= 0) = 1;
    skels{i}.nodesNumDataAll{1}(:, 3:5) = skels{i}.nodes{1}(:, 1:3);
    skels{i} = Skeleton.setParams4Dataset(skels{i}, 'ex145_ROI2017');
end
end
