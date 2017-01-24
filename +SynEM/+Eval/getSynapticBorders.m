function [ isSynaptic, group, unlabeledGT ] = getSynapticBorders( ...
    borders, segGT, overlapT, minOverlap )
%GETSYNAPTICBORDERS Determine the synaptic borders by labeling all borders
%as synaptic if the overlap with the ground truth segmentation.
% INPUT borders: [Nx1] border struct.
%       segGT: Logical or label array containing the ground truth
%           segmentation.
%       overlapT: (Optional) float
%           Fraction of a border that needs to be synaptic in order to label
%           the border as synaptic.
%           (Default: 0)
%       minOverlap: (Optional) double
%           Absolute number of minimal synaptic voxels in an interface to
%           be considered synaptic. If both overlapT and minOverlap are
%           specified than interfaces have to satisfy both requirements.
%           (Default: 0).
% OUTPUT isSynaptic: [Nx1] logical
%           Logical indices of the synaptic borders in borders.
%        group: [Nx1] int
%           Integer array grouping all borders belonging to the same ground
%           truth synapse.
%        unlabeledGT: [Nx1] int
%           IDs of ground truth clusters that did not label any interface.
%           If segGT is a logical matrix then the ids are w.r.t.
%           bwabeln(segGT)
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('overlapT','var') || isempty(overlapT)
    overlapT = 0;
end
if ~exist('minOverlap','var') || isempty(minOverlap)
    minOverlap = 0;
end

isSynaptic = false(length(borders),1);
group = zeros(length(borders),1,'uint16');

%convert segGT to label matrix if necessary
if islogical(segGT)
    L = bwlabeln(segGT);
else
    L = segGT;
    segGT = L > 0;
end

for i = 1:length(borders)
    isSynaptic(i) = (sum(segGT(borders(i).PixelIdxList)) > ...
        max(overlapT*length(borders(i).PixelIdxList),minOverlap));
    if isSynaptic(i)
        group(i) = firstEl(setdiff(unique(L(borders(i).PixelIdxList)),0));
    end
end
unlabeledGT = setdiff(setdiff(L(:),group(:)), 0);
end

function a = firstEl(A)
a = A(1);
end

