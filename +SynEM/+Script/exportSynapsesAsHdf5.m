% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
% See https://gitlab.mpcdf.mpg.de/connectomics/Benedikt/blob/b9b51fc9ef4261543f43e648a1b36661c1741923/+SynEM/+Training/+runConfigShaft/cnn17_2_synem_trShaft.m
inDir = '/gaba/u/bstaffle/data/SynEM/SynapseDetectionTrainingData_2_typeLabels';
outDir = '/home/amotta/Desktop';

info = Util.runInfo();
Util.showRunInfo(info);

%% Convert MAT to HDF5 files
clear cur*;

matFiles = dir(fullfile(inDir, '*.mat'));
matFiles(cat(1, matFiles.isdir)) = [];
matFiles = {matFiles.name};

%%
for curIdx = 1%:numel(matFiles)
    curMatFile = matFiles{curIdx};
    
    curInFilePath = fullfile(inDir, curMatFile);
    curIn = load(curInFilePath);
    
    curOut = struct;
    curOut.em = curIn.data.raw;
    curOut.segmentation = curIn.data.segments;
    
    curOut.interfaces = cellfun( ...
        @(voxelIds) voxelIds(:), ...
        curIn.data.interfaceSurfaceList(:), ...
        'UniformOutput', false);
    curOut.edges = curIn.data.neighborIDs;
    
    curOut.labels = curIn.interfaceLabels;
    curOut.labels(curOut.labels == 3) = 0;
    curOut.labels = categorical( ...
        curOut.labels, uint8(0:2), ...
        {'NoSynapse', 'PreToPostSynapse', 'PostToPreSynapse'});
    
    curOut.targetTypes = zeros(size(curOut.interfaces), 'uint8');
    curOut.targetTypes(curIn.typeLabels.idx( ...
        curIn.typeLabels.type == 'spine_head')) = 1;
    curOut.targetTypes(curIn.typeLabels.idx( ...
        curIn.typeLabels.type == 'spine_neck')) = 2;
    curOut.targetTypes(curIn.typeLabels.idx( ...
        curIn.typeLabels.type == 'shaft')) = 3;
    curOut.targetTypes(curIn.typeLabels.idx( ...
        curIn.typeLabels.type == 'soma')) = 4;
    curOut.targetTypes = categorical( ...
        curOut.targetTypes, uint8(0:4), ...
        {'NoSynapse', 'SpineHead', 'SpineNeck', 'Shaft', 'Soma'});
    
   [~, curOutFile] = fileparts(curMatFile);
    curOutFile = sprintf('%s.hdf5', curOutFile);
    curOutFile = fullfile(outDir, curOutFile);
    structToHdf5(curOutFile, '/', curOut);
    infoToHdf5(curOutFile, info);
end
