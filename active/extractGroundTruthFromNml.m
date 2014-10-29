function [labeledIdx, labels] = extractGroundTruthFromNml(seg, edges, skelParam)

skel = parseNml(skelParam.file);
skel = removeEmptySkeletons(skel);
skel = switchToLocalCoords_v2(skel, skelParam.bbox(:,1)');
skel = correctSkeletonsToBBox_v2(skel, skelParam.bbox(:,2)' - skelParam.bbox(:,1)');
% Find labels
equivMatrix = evalSegError(seg, skel, 1); 
[~, idx] = find(sum(equivMatrix,1) ~= 0);
labeledIdx = false(1,size(edges,1));
labels = [];
for i=1:size(edges,1)
	if any(edges(i,1) == idx) && any(edges(i,2) == idx)
		labeledIdx(i) = true;
		vec1 = equivMatrix(:,edges(i,1));
		vec2 = equivMatrix(:,edges(i,2));
		if any(vec1 & vec2)
			labels(end+1) = 1;
		else
			labels(end+1) = -1;
		end
	end
end
% Transpose in order to have same dimension as supervoxel weights
labeledIdx = labeledIdx';
labels = labels';

end

