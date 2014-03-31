function [ result, longSeg] = correspondenceFinder(p)
% correspondenceFinder: Generate struct with equivalence relations between
% segmentation cubes using 2 planes of overlap

% This assumes all bounding boxes of local segmentations have the same size:
cubeSize = p.local(1,1,1).bboxBig(:,2) - p.local(1,1,1).bboxBig(:,1);

% list of tuples of adjacent cubes
adjCubes = getOverlaps(size(p.local));

% For every tuple of adjacent cubes calculated correspondences
result = struct();
for i = 1:size(adjCubes,1)
    functionH{i} = @calculateLocalCorrepondences;
    inputCell{i} = {adjCubes(i,1:3), adjCubes(i,4:6), p.local(1,1,1).
    result(i).data  = globalize(adjCubes(i,1:3) ,adjCubes(i,4:6), param.cubeSize, param.Overlap, param.segPath);
end

end

%%
% Param File

% param.cubeLims = [ 1 4; 1 4; 1 4]
% param.cubeSize = [ 1024; 1024; 512]
% param.Overlap = 1
% param.segPath = 'I:\CortexConnectomics\Manuel\results\pipeline\20140312T141921\'
