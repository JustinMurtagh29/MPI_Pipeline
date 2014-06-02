function [skel, accStruct] = accumulateOverCubeBorders(skel, aStruct, nameForTree)
% Pass structure with field voxel (in global coordinates) positions, will be merged if overlapping, name for tree to add to skel

nrElements = length(aStruct);
adjList = cell(nrElements,1);
% Create symetric adjaceny list of elements
for i=1:nrElements-1
    for j=i+1:nrElements
	a = aStruct(i).voxel;
        b = aStruct(j).voxel;
	if ~isempty(intersect(a, b, 'rows')) 
           adjList{i}(end+1) = j;
	   adjList{j}(end+1) = i;
	end
    end   
end
% Find CC
components = bfs(1:nrElements, adjList);
for i=1:length(components)
    accStruct(i).voxel = unique(vertcat(aStruct(components{i}).voxel), 'rows');
    skel = addTree(skel, accStruct(i).voxel, [nameForTree num2str(i, '%.2i')]);
end

end
