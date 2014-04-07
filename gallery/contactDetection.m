function [skel, merger] = contactDetection(skeletons, outputDir)
% Pass exactly two file names in a cell array
% Modified Benjamins gallery code (quite a bit)

skel1 = parseNml(skeletons{1});
skel2 = parseNml(skeletons{2});
skel = mergeTrees(skel1, skel2, skeletons{1}, skeletons{2});
% for each node from both skeletons, find cube in which it resides
cubeCoords =  [floor((skel{1}.nodes(:,1:3)-1)./128); floor((skel{2}.nodes(:,1:3)-1)./128)];
% all unique cubes traversed by this skeleton 
uniqueCubeCoords = unique(cubeCoords, 'rows');
% contact counter
contactNumber = 0;
mergerNumber = 0;
for i = 1:size(uniqueCubeCoords,1)
    % Bounding box of current unique cube including overlap (hardcoded for 64 border)
    lowerLimit = uniqueCubeCoords(i,:).*128+1-64;
    upperLimit = lowerLimit+255;
    % find all nodes from skeletons that lie in the same cube including overlap(!)
    nodeIdx1 = all(bsxfun(@le, skel{1}.nodes(:,1:3), upperLimit) & bsxfun(@ge, skel{1}.nodes(:,1:3), lowerLimit),2);
    nodeIdx2 = all(bsxfun(@le, skel{2}.nodes(:,1:3), upperLimit) & bsxfun(@ge, skel{2}.nodes(:,1:3), lowerLimit),2);
    % Continue only if nodes in cube (including overlap) in other skeleton is present 
    if any(nodeIdx1) && any(nodeIdx2)
	% all nodes for this cube from both skeletons (with overlap for 2, without for 1)
        skel1nodes = skel{1}.nodes(nodeIdx1,1:3);
        skel2nodes = skel{2}.nodes(nodeIdx2,1:3);
	skel1nodesLocal = bsxfun(@minus, skel1nodes, lowerLimit-1);
	skel2nodesLocal = bsxfun(@minus, skel2nodes, lowerLimit-1);
    	% read local segmentation from disk
	cube = readKnossosCube('/nfs/bmo/mberning/20140310backup/mag1/', '100527_k0563_seg', uniqueCubeCoords(i,:), 'uint16', '', 'raw', 256);
	% get all local segmentation IDs of the nodes in the cube
	segIds1 = cube(sub2ind(size(cube),skel1nodesLocal(:,1),skel1nodesLocal(:,2),skel1nodesLocal(:,3)));
	uniqueSegIds1 = unique(segIds1);
	uniqueSegIds1(uniqueSegIds1 == 0) = [];
	segIds2 = cube(sub2ind(size(cube),skel2nodesLocal(:,1),skel2nodesLocal(:,2),skel2nodesLocal(:,3)));
	uniqueSegIds2 = unique(segIds2);
	uniqueSegIds2(uniqueSegIds2 == 0) = [];
	% Detect merger: winner takes all according to number of nodes in segment/supervoxel
	mergeSegId = intersect(uniqueSegIds1, uniqueSegIds2); 
	if ~isempty(mergeSegId);
	    display([num2str(uniqueCubeCoords(i,:), '%.2i,%.2i,%.2i') ' - Merger detected, ID(s): ' num2str(mergeSegId)]);
	    for j=1:length(mergeSegId)
		mergerNumber = mergerNumber + 1;
		merger(mergerNumber).cubeId = uniqueCubeCoords(i,:);
		merger(mergerNumber).id = mergeSegId(j);
		temp = regionprops(cube == mergeSegId(j), 'PixelList');
		merger(mergerNumber).voxel = bsxfun(@plus,temp.PixelList(:,[2 1 3]),lowerLimit-1);
		merger(mergerNumber).nodeCounts(1) = sum(segIds1 == mergeSegId(j));
		merger(mergerNumber).nodeCounts(2) = sum(segIds2 == mergeSegId(j));
		% Winner takes all, nobody takes anything in case of tie
		if merger(mergerNumber).nodeCounts(1) > merger(mergerNumber).nodeCounts(2)
		    uniqueSegIds2(uniqueSegIds2 == mergeSegId(j)) = [];
		elseif merger(mergerNumber).nodeCounts(2) > merger(mergerNumber).nodeCounts(1)
		    uniqueSegIds1(uniqueSegIds1 == mergeSegId(j)) = [];
		else
		    uniqueSegIds1(uniqueSegIds1 == mergeSegId(j)) = [];
		    uniqueSegIds2(uniqueSegIds2 == mergeSegId(j)) = [];
		end
	    end	
	end
	% Find contact area in segmentation
	cube1 = ismember(cube,uniqueSegIds1);
	cube2 = ismember(cube,uniqueSegIds2);
	contactCube = imdilate(cube1, ones(3,3,3)) & imdilate(cube2, ones(3,3,3));
	contactCube = imclose(contactCube, ones(5,5,5));
	props = regionprops(contactCube, {'PixelList'});
	for j=1:length(props)
		contactNumber = contactNumber + 1;
		contact(contactNumber) = props(j);
		contact(contactNumber).PixelList = bsxfun(@plus,contact(contactNumber).PixelList(:,[2 1 3]),lowerLimit-1);
	end
    end
end

% accumulate overlapping contacts
adjMatrix = zeros(contactNumber, 'uint8');
for i=1:contactNumber-1
    for j=i+1:contactNumber
	adjMatrix(i,j) = ~isempty(intersect(contact(i).PixelList, contact(j).PixelList, 'rows'));
    end
end

% add accumulated contacts as new tree to skeleton struct each
for i=1:contactNumber
	skel = addTree(skel, ['Contact ID ' num2str(i, '%.2i') ' of size ' num2str(contact(i).Area)], contact(i).Centroid);
end

% accumulate merger
for i=1:length(merger)

end

end

