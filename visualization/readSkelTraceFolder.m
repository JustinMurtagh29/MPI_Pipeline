function output = readSkelTraceFolder( path )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
output = struct('trace', {});
folderList = dir([path 'seed*']);
for i = 1:length(folderList)
    if isdir([path folderList(i).name])
        skelList = dir([path folderList(i).name '\*.nml']);
        for j = 1:length(skelList)
                output(i,j).trace = readKnossosNml([path folderList(i).name '\' skelList(j).name], 1);
        end
    end
end
end

