function job = correspondenceFinder(p)
% correspondenceFinder: Generate struct with equivalence relations between
% segmentation cubes using 2 planes of overlap

% list of tuples of adjacent cubes
adjCubes = getOverlaps(p.tiles);
inputCell = cell(size(adjCubes, 1), 1);

% For every tuple of adjacent cubes: calculate correspondences of segmentation IDs between cubes using overlap
for i = 1:size(adjCubes,1)
    inputCell{i} = {adjCubes(i,1:3), adjCubes(i,4:6), ...
        p.local(adjCubes(i,1),adjCubes(i,2),adjCubes(i,3)).tempSegFile, ...
        p.local(adjCubes(i,4),adjCubes(i,5),adjCubes(i,6)).tempSegFile, ...
        p.local(adjCubes(i,1),adjCubes(i,2),adjCubes(i,3)).bboxSmall, ...
        p.local(adjCubes(i,4),adjCubes(i,5),adjCubes(i,6)).bboxSmall, ...
        p.local(adjCubes(i,1),adjCubes(i,2),adjCubes(i,3)).bboxBig, ...
        p.local(adjCubes(i,4),adjCubes(i,5),adjCubes(i,6)).bboxBig, ... 
        p.correspondence.overlap, p.correspondence.saveFolder};
end

functionH = @calculateLocalCorrespondences;
job = startCPU(functionH, inputCell, 'correspondence', 12, 100);

end

