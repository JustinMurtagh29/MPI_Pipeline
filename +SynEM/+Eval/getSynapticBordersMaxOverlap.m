function isSynaptic = getSynapticBordersMaxOverlap(borders, segGT)
%GETSYNAPTICBORDERSMAXOVERLAP Label the borders with maximal overlap for
% each segmentation connected components.
% INPUT borders: [Nx1] border struct.
%       segGT: Logical or label array containing the ground truth
%           segmentation.
% OUTPUT isSynaptic: [Nx1] logical
%           Logical vector indicating for each border whether it is
%           synpatic. For each segmentation connected component only the
%           border with maximal overlap (in total number of overlapping
%           voxels) is labeled.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if islogical(segGT)
    L = bwlabeln(segGT);
else
    L = segGT;
end
synIds = setdiff(L(:), 0);
numSyn = length(synIds);

borderSegIds = cellfun(@(x)L(x), {borders.PixelIdxList}, 'uni', 0);
isSynaptic = false(length(borders), 1);
for i = 1:numSyn
    [~,idx] = max(cellfun(@(x)sum(x == synIds(i)), borderSegIds));
    isSynaptic(idx) = true;
end


end
