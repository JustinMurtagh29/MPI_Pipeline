% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/tmpscratch/amotta/l23/2018-10-09-mrnet-pipeline-run';
graphFile = '/gaba/u/amotta/l23/2018-10-11-hierarchical-agglomeration/graph.mat';

scoreThresh = -0.58;
distThreshNm = 10E3;

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
param.experimentName = '2018-10-10_ex144_st08x2_mrnet';

maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);

graphT = load(graphFile);
graphT = graphT.graph;

%% Limit graph (mainly to reduce RAM usage)
graphT = graphT(graphT.score > scoreThresh, :);

%% Build nuclear agglomerates
curSeedPos = [3669, 4737, 3983] + 1;
curSeedSegId = Seg.Global.getSegIds(param, curSeedPos);

curBox = round(distThreshNm ./ param.raw.voxelSize(:));
curBox = transpose(curSeedPos(:) + [-1, +1] .* curBox);

curCandSegIds = find( ...
    all(segPoints >= curBox(1, :), 2) ...
  & all(segPoints <= curBox(2, :), 2));
assert(any(curCandSegIds == curSeedSegId));

curGraphT = graphT;
[curMask, curEdgeIds] = ismember(graphT.edge, curCandSegIds);
curGraphT = curGraphT(all(curMask, 2), :);
curGraphT.edge = curEdgeIds(all(curMask, 2), :);

curGraph = graph( ...
    curGraphT.edge(:, 1), ...
    curGraphT.edge(:, 2), ...
    curGraphT.score, ...
    numel(curCandSegIds));
curComps = conncomp(curGraph);

curAggloId = curComps(curCandSegIds == curSeedSegId);
curAgglo = curCandSegIds(curComps == curAggloId);

%% Inspect
curSkel = skeleton();
curSkel = Skeleton.setParams4Pipeline(curSkel, param);
curSkel = curSkel.addTree('Nucleus', segPoints(curAgglo, :));
