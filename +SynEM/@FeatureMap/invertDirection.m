function X = invertDirection(obj, X)
%INVERTDIRECTION Invert direction of interface feature matrix.
% INPUT X: [NxM] float where M = obj.numFeaturesSelected
% OUTPUT X: [NXM] float input feature matrix with columns interchanged that
%           belong to different subsegments.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if isempty(obj.selectedFeat)
    error('Run setSelectedFeat first.')
end

numFeatText = length(obj.featTexture);
invIdx = cell(numFeatText + length(obj.featShape),1);
count = 0;
for i = 1:numFeatText
    % calculate inversion index for texture features
    invIdx{i} = double(obj.selectedFeat{i});
    numFeat = sum(invIdx{i}(:));
    invIdx{i}(obj.selectedFeat{i}) = (1:numFeat) + count;
    invIdx{i} = invIdx{i}(:,[1 reshape([3:2:2*(obj.numSubvolumes+1); ...
            2:2:2*(obj.numSubvolumes)],1,[])],:);
    invIdx{i} = invIdx{i}(obj.selectedFeat{i});
    count = count + numFeat;
end

for i = 1:length(obj.featShape)
    %get direction dependency of current features
    tmp = obj.featShapeDirDep{i}(obj.selectedFeat{i + numFeatText});

    %renumber indices starting from one (required if not all features are
    %selected)
    [~,~,tmp] = unique(tmp);
    invIdx{i + numFeatText} = tmp + count;
    count = count + length(tmp);
end
invIdx = cell2mat(invIdx);
X = X(:,invIdx);
end
