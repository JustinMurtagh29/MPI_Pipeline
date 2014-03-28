function [ result, longSeg] = getData()
%UNTITLED Generate struct with equivalence relations between
%segmentation overlaps

param = load('P:\globalization\param.mat'); %for CubeLims, CubeSize, Overlap, segmentation_path
param = param.param;

%list of tuples of adjacent cubes
adjCubes = getOverlaps(param.cubeLims);

%Load cube-segmentations and create field in output struct for every tuple
%of adjacent cubes

result = struct();
for i = 1:1%size(adjCubes,1)
    result(i).data  = globalize(adjCubes(i,1:3) ,adjCubes(i,4:6), param.cubeSize, param.Overlap, param.segPath);
end


%%
% Param File

% param.cubeLims = [ 1 4; 1 4; 1 4]
% param.cubeSize = [ 1024; 1024; 512]
% param.Overlap = 1
% param.segPath = 'I:\CortexConnectomics\Manuel\results\pipeline\20140312T141921\'