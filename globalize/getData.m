function [ segmentOverlap ] = getData( cubeLims )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%list of tuples of adjacent cubes
adjCubes = getOverlaps(cubeLims);

segmentOverlap = struct();
for i = 1:size(adjCubes,1)
    % data : idObject1  idObject2  voxelSize1  voxelSize2  NrOfCommonVoxels
    segmentOverlap(i).data = globalize(adjCubes(i,1:3),adjCubes(i,4:6));
    segmentOverlap(i).cubeCoords1 = adjCubes(i,1:3);
    segmentOverlap(i).cubeCoords2 = adjCubes(i,4:6);
end