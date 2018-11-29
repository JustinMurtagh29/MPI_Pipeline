function job = correspondenceFinderGlobal(p)
% correspondenceFinder: Generate struct with equivalence relations between
% segmentation cubes using adjacent planes across cube border

% Save to new folder for now, should be integrated into pipeline at some point
saveFolder = p.correspondence.saveFolder;

if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

% list of tuples of adjacent cubes
adjCubes = getOverlaps(p.tiles);

% For every tuple of adjacent cubes: calculate correspondences of segmentation IDs between cube
inputCell = cell(size(adjCubes, 1), 1);
for i = 1:size(adjCubes, 1)
    curIdsOne = adjCubes(i, 1:3);
    curIdsTwo = adjCubes(i, 4:6);
    
    inputCell{i} = { ...
        p.seg, curIdsOne, curIdsTwo, ...
        p.local(curIdsOne(1), curIdsOne(2), curIdsOne(3)).bboxSmall, ...
        p.local(curIdsTwo(1), curIdsTwo(2), curIdsTwo(3)).bboxSmall, ...
        saveFolder};
end

functionH = @calculateGlobalCorrespondences;
job = startCPU(functionH, inputCell, 'correspondence', 12, 50, 0);

end

