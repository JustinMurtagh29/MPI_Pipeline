% This script calculates the isosurfaces of whole cells at three different
% stages of spine attachment:
%
% * before any spine attachment
% * after automated spine attachment
% * after automated + manual spine attachment
%
% At a later stage, it might also be useful for this script to export the
% skeleton tracings for manually attached spines.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outputDir = '/tmpscratch/amotta/l4/2018-05-10-whole-cell-isosurfaces-spine-evolution';

dendFiles = struct;
% Prior to spine attachment
dendFiles(1).file = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3.mat');
dendFiles(1).tag = 'pre';
% After automated spine attachment
dendFiles(2).file = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');
dendFiles(2).tag = 'auto';
% After automated + manual spine attachment
dendFiles(3).file = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto-and-manual.mat');
dendFiles(4).tag = 'full';

% Segmentation
segParam = struct;
segParam.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
segParam.backend = 'wkwrap';

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;
param.seg = segParam;

%% Processing dendrite stages
for curFile = dendFiles
    Util.log('Processing "%s" state (%s)', curFile.tag, curFile.file);
    curDir = fullfile(outputDir, curFile.tag);
    curDend = load(curFile.file);
    
    curWcT = table;
    curWcT.id = curDend.idxWholeCells;
    curWcT.agglo = curDend.dendrites;
    
    % Only keep whole cells
    curWcT(~curWcT.id, :) = [];
    curWcT = sortrows(curWcT, 'id');
    
    % Convert to segment equivalence classes
    curWcT.agglo = Agglo.fromSuperAgglo(curWcT.agglo, true);
    
    % Save meta data
    curOut = struct;
    curOut.file = curFile;
    curOut.wholeCells = curWcT;
    curOut.info = info;
    
    %% Writing results
    mkdir(curDir);
    
    curOutFile = fullfile(curDir, 'info.mat');
    Util.saveStruct(curOutFile, curOut);
    
    Visualization.exportAggloToAmira( ...
        param, curWcT.agglo, curDir, ...
        'smoothSizeHalf', 4, ...
        'smoothWidth', 8, ...
        'reduce', 0.05);
end

%% Done
Util.log('Done');
