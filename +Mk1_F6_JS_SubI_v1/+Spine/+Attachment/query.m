% Build queries for manual spine head attachment.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
configName = 'ex144_08x2_mrNet';
dendFile = '/tmpscratch/amotta/l23/2019-11-18-axon-and-dendrite-agglomerates/20191120T024253_results_auto-spines_v2.mat';
outFile = '/tmpscratch/amotta/l23/2020-01-14-manual-spine-head-attachment/task-definition_v1.mat';

% NOTE(amotta): If set, an NML file with a random subset of unattached
% spine heads will be written to this location.
skelFile = '';

% As in L4 paper
marginNm = 3000;

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
config = loadConfig(configName);
param = config.param;

voxelSize = param.raw.voxelSize;

segSizes = Seg.Global.getSegToSizeMap(param);
segCentroids = Seg.Global.getSegToCentroidMap(param);
segPoints = Seg.Global.getSegToPointMap(param);

[shAgglos, attached] = Util.load(dendFile, 'shAgglos', 'attached');

%% Calculate spine head positions
shPos = nan(numel(shAgglos), 3);
wmean = @(v, w) sum(v .* w(:), 1) / sum(w(:));

for curIdx = 1:numel(shAgglos)
    curSegIds = shAgglos{curIdx};
    
    curCentroid = wmean( ...
        segCentroids(curSegIds, :), ...
        segSizes(curSegIds));
    
   [~, curRelSegId] = pdist2( ...
       voxelSize .* segPoints(curSegIds, :), ...
       voxelSize .* curCentroid, ...
       'squaredeuclidean', 'Smallest', 1);
    
    curShPos = curSegIds(curRelSegId);
    curShPos = segPoints(curShPos, :);
    
    shPos(curIdx, :) = curShPos;
end

%% Inspect random subset in webKnossos
clear cur*;
rng(0);

if ~isempty(skelFile)
    curIds = find(not(attached));
    curRandIds = randperm(numel(curIds));
    curIds = curIds(curRandIds(1:25));

    curNumDigits = ceil(log10(1 + numel(curIds)));

    skel = skeleton();
    skel = Skeleton.setParams4Pipeline(skel, param);
    skel = Skeleton.setDescriptionFromRunInfo(skel, info);

    for curIdx = 1:numel(curIds)
        curId = curIds(curIdx);
        curPos = shPos(curId, :);

        curName = sprintf( ...
            '%0*d. Spine head %d', ...
            curNumDigits, curIdx, curId);

        skel = skel.addTree(curName, curPos);
    end

    skel.write(skelFile);
end

%% Generate queries
clear cur*;

curMarginVx = ceil(marginNm ./ voxelSize);

curCoreMask = ...
    all((shPos - param.bbox(:, 1)') > curMarginVx, 2) ...
  & all((param.bbox(:, 2)' - shPos) > curMarginVx, 2);

unattached = not(attached);
unattachedCore = unattached & curCoreMask;

%% Save results
clear cur*;

curQueryPos = shPos(unattached, :);

out = struct;
out.info = info;
out.shQueryPos = curQueryPos;
Util.saveStruct(outFile, out);
Util.protect(outFile);
