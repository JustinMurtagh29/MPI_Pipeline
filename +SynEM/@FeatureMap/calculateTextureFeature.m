function [ X ] = calculateTextureFeature( obj, vols, raw, i, ignoreBorder )
%CALCULATETEXTUREFEATURE Wrapper for texture feature calculation.
% INPUT vols, raw, ignoreBorder: see obj.calculate
%       i: Number of the texture feature to calculate.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%calculate the feature
feat = obj.featTexture{i}.calc(raw);

%convert to cell if only one channel (for loop below)
if ~iscell(feat)
    feat = {feat};
end

%crop feature to valid bounding box
if ~ignoreBorder
    b = obj.border./2;
    feat = cellfun(@(x)x(b(1) + 1:end - b(1), b(2) + 1: end - b(2), ...
        b(3) + 1: end - b(3)),feat,'UniformOutput',false);
end

%preallocate output
X = zeros(size(vols,1),sum(obj.selectedFeat{i}(:)),'like', feat{1});

%calculate row indices for subvolumes
rowIdx = sum(obj.selectedFeat{i},1);
rowIdx = [0; cumsum(rowIdx(:))];
nV = 2*obj.numSubvolumes + 1;

%pool from filter reponse
for j = 1:obj.featTexture{i}.numChannels
    for k = 1:nV %loop over subvolumes
        X(:,rowIdx((j-1)*nV + k) + 1:rowIdx((j-1)*nV + k + 1)) = ...
            cell2mat( cellfun(@(x)obj.calcSumStatsInternal ...
            (feat{j}(x),obj.selectedFeat{i}(:,k,j)), ...
            vols(:,k),'UniformOutput',false));
    end
end
end
