% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
clear;

minCv = 0.5;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');
asiRunId = '20190227T082543';

outDir = '/home/amotta/Desktop';

info = Util.runInfo();
Util.showRunInfo(info);

%% Load axon-spine interfaces
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, ~, connFile] = connectEM.Consistency.loadConnectomePaper(param);

shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

axons = load(conn.info.param.axonFile);
axons = axons.axons(axons.indBigAxons);
dendrites = load(conn.info.param.dendriteFile);
dendrites = dendrites.dendrites;

[curDir, curAsiFile] = fileparts(connFile);
curAsiFile = sprintf('%s__%s_asiT.mat', curAsiFile, asiRunId);
curAsiFile = fullfile(curDir, curAsiFile);

asiT = load(curAsiFile);
asiT = asiT.asiT;

segPoints = Seg.Global.getSegToPointMap(param);

%% Building pairs
clear cur*;

asiT = asiT(asiT.area > 0, :);
asiT = connectEM.Consistency.Calibration.apply(asiT);

asiT = asiT(asiT.type == 'PrimarySpine', :);
curExcAxonClasses = {'Corticocortical', 'Thalamocortical'};
asiT = asiT(ismember(asiT.axonClass, curExcAxonClasses), :);

curPairs = struct('synIds', 1:height(asiT));
curPairs = connectEM.Consistency.buildPairConfigs(asiT, curPairs);
curPairs = curPairs(1);

sasdT = table;
sasdT.asiIds = curPairs.synIdPairs;
sasdT.synIds = asiT.id(sasdT.asiIds);
sasdT.shIds = asiT.shId(sasdT.asiIds);
sasdT.areas = asiT.area(sasdT.asiIds);
sasdT.cv = std(sasdT.areas, 0, 2) ./ mean(sasdT.areas, 2);
sasdT.axonId = asiT.preAggloId(sasdT.asiIds(:, 1));
sasdT.dendId = asiT.postAggloId(sasdT.asiIds(:, 1));

%% Export NMLs for high-CV pairs
clear cur*;
rng(0);

curSasdIds = find(sasdT.cv > minCv);
curSasdIds = curSasdIds(randperm(numel(curSasdIds)));
curSasdIds = curSasdIds(1:20);

curNumDigits = ceil(log10(1 + numel(curSasdIds)));

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = Skeleton.setDescriptionFromRunInfo(skel, info);

for curSasdIdx = 1:numel(curSasdIds)
    curSasdId = curSasdIds(curSasdIdx);
    curSasdT = sasdT(curSasdId, :);
    
    curAxon = axons(conn.axonMeta.parentId(curSasdT.axonId));
    curDend = dendrites(conn.denMeta.parentId(curSasdT.dendId));
    
    curShAgglos = shAgglos(curSasdT.shIds);
    curShSegIds = cellfun(@(segIds) segIds(1), curShAgglos);
   [~, curShNodeIds] = ismember(curShSegIds, curDend.nodes(:, 4));
   
    curComments = repelem( ...
        {''}, numel(curDend.nodes(:, 1)), 1);
    curComments(curShNodeIds) = arrayfun( ...
        @(synId) sprintf('Synapse %d', synId), ...
        curSasdT.synIds, 'UniformOutput', false);
    
    curSkel = skel;
    curSkel = curSkel.addTree( ...
        sprintf('Axon %d', curSasdT.axonId), ...
        curAxon.nodes(:, 1:3), curAxon.edges);
    curSkel = curSkel.addTree( ...
        sprintf('Dendrite %d', curSasdT.dendId), ...
        curDend.nodes(:, 1:3), curDend.edges, ...
        [], [], curComments);
    
    curNmlName = sprintf( ...
        '%0*d_axon-%d_dendrite_%d.nml', curNumDigits, ...
        curSasdIdx, curSasdT.axonId, curSasdT.dendId);
    curSkel.write(fullfile(outDir, curNmlName));
end
