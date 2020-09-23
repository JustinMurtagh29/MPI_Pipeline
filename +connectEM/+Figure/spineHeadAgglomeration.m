% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
shAggloFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');
outDir = '/home/amotta/Desktop';

% See https://gitlab.mpcdf.mpg.de/connectomics/amotta/blob/b3d6b7be4f876c54e4261026d5ab47d4edd49d89/matlab/+L4/+Figure/highResEmSamples.m
% Range from Benedikt, used for SynEM paper.
emRange = [60, 180];

configs = struct;
configs(1).shId = 1626;
configs(1).pos = [2957, 3384, 349] + 1;
configs(1).fovUm = 1.2;

configs(2).shId = 213960;
configs(2).pos = [4475, 2021, 1497] + 1;
configs(2).fovUm = 1.2;

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

graphT = Graph.load(param.saveFolder);

shAgglos = load(shAggloFile);
shAgglos = shAgglos.shAgglos;

maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);

%% Quantitative evaluation
clear cur*;

% Fraction of spine head with multiple segments
curMultiFrac = cellfun(@numel, shAgglos);
curMultiFrac = mean(curMultiFrac > 1) %#ok

% Approx. fraction of spine heads that were split by cube boundaries
% NOTE(amotta): This is only an approximation - albeit a very good one -
% because we are also counting low-probability edges that might not have
% been used during spine head agglomeration.
%   A quick check showed that more than 99.7 % of edges collected this have
% actually been used during agglomeration.
curShLUT = Agglo.buildLUT(maxSegId, shAgglos);
graphT.edges = curShLUT(graphT.edges);

graphT = graphT(all(graphT.edges, 2), :);
graphT = graphT(graphT.edges(:, 1) == graphT.edges(:, 2), :);

curIntraCubeSplitMask = accumarray( ...
    graphT.edges(:, 1), graphT.borderIdx, ...
   [numel(shAgglos), 1], @any, false);
curIntraCubeSplitFrac = mean(curIntraCubeSplitMask) %#ok

%% Export single- and multi-segment spine head examples
clear cur*;
rng(0);

curSample = @(a) feval( ...
    @(b) b(1:min(20, numel(b))), ...
    a(randperm(numel(a))));

curSegCount = cellfun(@numel, shAgglos);
curSingleIds = curSample(find(curSegCount == 1)); %#ok
curMultiIds = curSample(find(curSegCount > 1)); %#ok

curIds = cat(1, curSingleIds(:), curMultiIds(:));
curAgglos = shAgglos(curIds);

curPoints = cellfun( ...
    @(ids) segPoints(ids, :), ...
    curAgglos, 'UniformOutput', false);

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = Skeleton.setDescriptionFromRunInfo(skel, info);
skel = Skeleton.fromMST(curPoints, param.raw.voxelSize, skel);

curNames = arrayfun( ...
    @(idx, id) sprintf( ...
        '%0*d. Spine head %d', ...
        ceil(log10(1 + numel(curIds))), idx, id), ...
    reshape(1:numel(curIds), [], 1), curIds(:), ...
    'UniformOutput', false);
skel.names = curNames;

%% Illustrate spine head agglomeration
clear cur*;

curColors = get(groot, 'defaultAxesColorOrder');

for curConfig = configs
    curBoxSize = 1E3 * curConfig.fovUm;
    curBoxSize = round(curBoxSize ./ param.raw.voxelSize);
    curBoxSize(end) = 1;
    
    curBox = curConfig.pos;
    curBox = round(curBox(:) - curBoxSize(:) / 2);
    curBox = curBox + [zeros(3, 1), curBoxSize(:) - 1];
    
    curRaw = transpose(loadRawData(param.raw, curBox));
    curSeg = transpose(loadSegDataGlobal(param.seg, curBox));
    
    curRaw = (double(curRaw) - emRange(1)) / diff(emRange);
    curRaw = min(max(curRaw, 0), 1);
    curRaw = repmat(curRaw, 1, 1, 3);
    
    curSegIds = shAgglos{curConfig.shId};
   [curMask, curSeg] = ismember(curSeg, curSegIds);
    
    curRaw = reshape(curRaw, [], 3);
    curRaw(curMask, :) = ...
        0.65 * curRaw(curMask, :) ...
      + 0.35 .* curColors(curSeg(curMask), :);
    curRaw = reshape(curRaw, [curBoxSize(1:2), 3]);
    
    curFig = figure();
    curAx = axes(curFig);
    imshow(curRaw, 'Parent', curAx);
    
    curFig.Position(3:4) = 300;
    curAx.Position = [0, 0, 1, 1];
    
    curOutName = sprintf('sh-agglo-%d.png', curConfig.shId);
    imwrite(curRaw, fullfile(outDir, curOutName));
end
