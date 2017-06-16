function job = correspondenceFinderGlobal(p)
% correspondenceFinder: Generate struct with equivalence relations between
% segmentation cubes using adjacent planes across cube border

% Save to new folder for now, should be integrated into pipeline at some point
saveFolder = [p.saveFolder 'correspondencesNew' filesep];
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

% list of tuples of adjacent cubes
adjCubes = getOverlaps(p.tiles);

% For every tuple of adjacent cubes: calculate correspondences of segmentation IDs between cube
for i = 1:size(adjCubes,1)
    inputCell{i} = {p.seg, adjCubes(i,1:3), adjCubes(i,4:6), ...
        [p.local(adjCubes(i,1),adjCubes(i,2),adjCubes(i,3)).saveFolder 'segGlobal.mat'], ...
        [p.local(adjCubes(i,4),adjCubes(i,5),adjCubes(i,6)).saveFolder 'segGlobal.mat'], ...
        p.local(adjCubes(i,1),adjCubes(i,2),adjCubes(i,3)).bboxSmall, ...
        p.local(adjCubes(i,4),adjCubes(i,5),adjCubes(i,6)).bboxSmall, ...
        saveFolder};
end
functionH = @calculateGlobalCorrespondences;
job = startCPU(functionH, inputCell, 'correspondence', 12, 100, 0);

end

