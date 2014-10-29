function [edgesNew, border, weights] = borderCalculation(parameter, edges)

% Load data
seg = loadSegData(parameter.seg.root, parameter.seg.prefix, parameter.bboxBig);
% Preallocate double size to avoid adding on the fly
edgesNew = zeros(size(edges,1)*2,2);
weights = zeros(size(edges,1)*2,5);
% Create border structure
border = struct('Area', {}, 'Centroid', {}, 'PixelIdxList', {});
% Find border(s) for each edge
idx = 1;
for i=1:size(edges,1)
	% Calculate border between objects
	object = seg == edges(i,1);
	otherObject = seg == edges(i,2);
	borderArray = imdilate(object, ones(3,3,3)) & imdilate(otherObject, ones(3,3,3));
	borderProp = regionprops(borderArray, {'Area' 'Centroid' 'PixelIdxList'});
	for j=1:length(borderProp)
		border(idx) = borderProp(j);
		edgesNew(idx,:) = edges(i,:);
		%% Segmentation based feature calculation (more needed?)
		properties = {'Centroid' 'PixelList', 'Area'};
		propObject = regionprops(object, properties);
		propOtherObject = regionprops(otherObject, properties);
		pc = princomp(propObject.PixelList, 'econ');
		directionObject = pc(:,1);
		pc = princomp(propOtherObject.PixelList, 'econ');
		directionOtherObject = pc(:,1);
		% Feature: border pixels between objects
		weights(idx,1) = length(borderProp(j).PixelIdxList);
		% Feature: size of objects
		weights(idx,2) = min(propObject.Area,propOtherObject.Area);
		weights(idx,3) = max(propObject.Area,propOtherObject.Area);
		% Feature: scalar product between major axes of objects
		weights(idx,4) = directionObject' * directionOtherObject;
		% Feature: distance between CoM of objects
		weights(idx,5) = norm(propObject.Centroid-propOtherObject.Centroid);
		idx = idx + 1;
	end
end

% Remove unnecessary preallocation
edgesNew(idx:end,:) = [];
weights(idx:end,:) = [];

end

