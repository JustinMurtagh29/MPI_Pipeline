% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/tmpscratch/amotta/l23/2018-10-09-mrnet-pipeline-run';
graphFile = '/gaba/u/amotta/l23/2018-10-11-hierarchical-agglomeration/graph.mat';

seedNmlFile = fileparts(mfilename('fullpath'));
seedNmlFile = fullfile(seedNmlFile, '+Data', 'nuclei.nml');

scoreThresh = 0.5;
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

%% Load seed positions and segments
clear cur*;

curNml = slurpNml(seedNmlFile);
curNodes = NML.buildNodeTable(curNml);
curNodes.coord = curNodes.coord + 1;
curTrees = NML.buildTreeTable(curNml);

% Santity checks
assert(isequal(sort(curTrees.id), sort(curNodes.treeId)));
assert(isequal(sort(curTrees.id), unique(curTrees.id)));

nucleusT = table;
nucleusT.label = categorical(lower(curTrees.name));
[~, nucleusT.nodePos] = ismember(curTrees.id, curNodes.treeId);
nucleusT.nodePos = curNodes.coord(nucleusT.nodePos, :);

nucleusT.seedSegId = ...
    Skeleton.getSegmentIdsOfNodes( ...
        param, nucleusT.nodePos, 26);
nucleusT.seedSegId = num2cell(nucleusT.seedSegId, 2);

nucleusT.seedSegId = cellfun( ...
    @(ids) mode(nonzeros(ids)), nucleusT.seedSegId);
assert(all(nucleusT.seedSegId > 0));

%% Build nuclear agglomerates
clear cur*;

nucleusT.agglo(:) = {[]};

for curId = 1:height(nucleusT)
    disp(curId);
    
    curSeedPos = nucleusT.nodePos(curId, :);
    curSeedSegId = nucleusT.seedSegId(curId);
    
    curBox = round(distThreshNm ./ param.raw.voxelSize(:));
    curBox = transpose(curSeedPos(:) + [-1, +1] .* curBox);

    curCandSegIds = find( ...
        all(segPoints >= curBox(1, :), 2) ...
      & all(segPoints <= curBox(2, :), 2));
    assert(any(curCandSegIds == curSeedSegId));

    curGraphT = graphT;
   [~, curEdgeIds] = ismember(graphT.edge, curCandSegIds);
   
    curGraphT.edge = curEdgeIds;
    curGraphT = curGraphT(all(curEdgeIds, 2), :);

    curGraph = graph( ...
        curGraphT.edge(:, 1), ...
        curGraphT.edge(:, 2), ...
        curGraphT.score, ...
        numel(curCandSegIds));
    curComps = conncomp(curGraph);

    curAggloId = curComps(curCandSegIds == curSeedSegId);
    curAgglo = curCandSegIds(curComps == curAggloId);
    nucleusT.agglo{curId} = curAgglo(:);
end

%% Inspect
clear cur*;

curSkel = skeleton();
curSkel = Skeleton.setParams4Pipeline(curSkel, param);

curDigits = ceil(log10(1 + height(nucleusT)));
curUp = @(s) [upper(s(1)), s(2:end)];

for curId = 1:height(nucleusT)
    curName = sprintf( ...
        'Nucleus %0*d. %s', curDigits, ...
        curId, curUp(char(nucleusT.label(curId))));
    curSkel = curSkel.addTree( ...
        curName, segPoints(nucleusT.agglo{curId}, :));
end
