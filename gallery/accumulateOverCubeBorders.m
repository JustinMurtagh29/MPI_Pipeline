function [skel, accStruct] = accumulateOverCubeBorders(skel, aStruct, nameForTree)
% Pass structure with field voxel (in global coordinates) positions, will be merged if overlapping, name for tree to add to skel

nrElements = length(aStruct);
adjList = cell(nrElements,1);
% Create symetric adjaceny list of elements
for i=1:nrElements-1
    for j=i+1:nrElements
	% Somehow this took longest to fix, MATLAB intersect with 'rows' only works if matrices have same number of colums, no fuck (warning/error) was given
	% This will make them equal length, but is pretty ugly, maybe if construct with padding not duplication better?
	a = zeros(max(size(aStruct(i).voxel,1),size(aStruct(j).voxel,1)),3);
	b = zeros(size(a));
	a(1:size(aStruct(i).voxel,1),:) = aStruct(i).voxel;
        b(1:size(aStruct(j).voxel,1),:) = aStruct(j).voxel;
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
