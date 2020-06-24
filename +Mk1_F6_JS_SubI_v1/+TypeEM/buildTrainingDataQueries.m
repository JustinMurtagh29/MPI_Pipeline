% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
datasetName = 'Mk1_F6_JS_SubI_v1';
%% Loading data
rootDir = ['/tmpscratch/sahilloo/data/' datasetName '/pipelineRun_mr2e_wsmrnet/'];
load(fullfile(rootDir,'allParameter.mat'))
param = p;

% output folder for saving new agglomeration state
outDir = fullfile(p.saveFolder, 'tracings', 'typeEM');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% NOTE(amotta): Random boxes #1, #14, and #5 are the ones handed out to
% HiWis for volume annotations. I've also used them to generate ConnectEM
% training data.
%   Subsequently, the remaining random boxes (as generated by
% L23.Scripts.buildRandomBoxes 7a79056c3c7dc3fb96e3f83314d12b7596a2fdda)
% were appended to this list. Some random boxes (e.g., number 4) were
% skipped if they were dominated by somatic volume.

boxesWk = [ ...
    9787, 10869, 2530, 267, 267, 108; ...
    10299, 12387, 2937, 267, 267, 108;...
    6744, 8960, 1544,267, 267, 108;...
    15456, 7492, 2245,267, 267, 108;...
    7316, 13403, 3211,267, 267, 108; ...
    9819, 13186, 1566,267, 267, 108; ...
    14785, 10056, 3347,267, 267, 108; ...
    12266, 12144, 1894,267, 267, 108; ...
    5020, 16781, 2255, 267, 267, 108; ...
    9357, 6399, 2176, 267, 267, 108;
    12983, 13638, 1771, 267, 267, 108; ...
    13475, 5814, 3377, 267, 267, 108; ,...
    14760, 14017, 3950, 267, 267, 108];

sampleRange = [1, 200];

info = Util.runInfo();
Util.showRunInfo(info);

segPoints = Seg.Global.getSegToPointMap(param);

%% Sample edges
clear cur*:

for curBoxIdx = 1:size(boxesWk, 1)
    curBoxWk = boxesWk(curBoxIdx, :);
    curBox = Util.convertWebknossosToMatlabBbox(curBoxWk);
    
    curSegIds = loadSegDataGlobal(param.seg, curBox);
    curSegIds = reshape(setdiff(curSegIds, 0), [], 1);
    
    curRange = min(sampleRange, numel(curSegIds));
    curRange = curRange(1):curRange(2);
    
    rng(0);
    curRandIds = randperm(numel(curSegIds));
    curRandIds = curSegIds(curRandIds(curRange));
    
    curDigits = ceil(log10(1 + numel(curRandIds)));
    
    curNames = arrayfun( ...
        @(idx, id) sprintf('%0*d. Segment %d', curDigits, idx, id), ...
        1:numel(curRandIds), reshape(curRandIds, 1, []), ...
        'UniformOutput', false);
    curNodes = segPoints(curRandIds, :);
    
    curSkel = skeleton();
    curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    curSkel = Skeleton.setDescriptionFromRunInfo(curSkel, info);
    
    curSkel = curSkel.addNodesAsTrees(curNodes);
    curSkel.names = reshape(curNames, [], 1);
    
    curSkel.write(fullfile(outDir, sprintf('box-%d.nml', curBoxIdx)));
end