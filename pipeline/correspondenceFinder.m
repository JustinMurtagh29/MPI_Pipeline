function [ result, longSeg] = correspondenceFinder(p)
% correspondenceFinder: Generate struct with equivalence relations between
% segmentation cubes using 2 planes of overlap

% list of tuples of adjacent cubes
adjCubes = getOverlaps(size(p.local));

% For every tuple of adjacent cubes: calculate correspondences of segmentation IDs between cubes using overlap
for i = 1:size(adjCubes,1)
    inputCell{i} = {adjCubes(i,1:3), adjCubes(i,4:6), p.local(adjCubes(i,1),adjCubes(i,2),adjCubes(i,3)).segFile, p.local(adjCubes(i,4),adjCubes(i,5),adjCubes(i,6)).segFile, ...
	p.tileBorder, p.correspondence.overlap, p.correspondence.saveFolder};
end

functionH = @calculateLocalCorrespondences;
startCPU(functionH, inputCell, 'correspondence');

end

