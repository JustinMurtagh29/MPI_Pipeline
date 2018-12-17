% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, connFile] = ...
    connectEM.Consistency.loadConnectome(param);

% Loading spine head agglomerates
shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

% Loading augmented graph
graph = Graph.load(rootDir);
graph(~graph.borderIdx, :) = [];

borderAreas = fullfile(rootDir, 'globalBorder.mat');
borderAreas = load(borderAreas, 'borderArea2');
borderAreas = borderAreas.borderArea2;

graph.borderArea = borderAreas(graph.borderIdx);
clear borderAreas;

%% Build axon-spine interface areas
asiT = ...
    connectEM.Connectome.buildAxonSpineInterfaces( ...
        param, graph, shAgglos, conn, syn);
asiT = asiT(asiT.type == 'PrimarySpine', :);
asiT.relId = reshape(1:height(asiT), [], 1);

%% Build connectivity matrix
clear cur*;

[pairT, ~, curIds] = unique(asiT(:, { ...
    'preAggloId', 'postAggloId'}), 'rows');
pairT.axonClass = conn.axonMeta.axonClass(pairT.preAggloId);
pairT.targetClass = conn.denMeta.targetClass(pairT.postAggloId);
pairT.asiIds = accumarray(curIds, asiT.relId, [], @(ids) {ids});

%% Actually analyse weight matrix
clear cur*;

curAxonClasses = {'Corticocortical', 'Thalamocortical'};
curTargetClasses = {};

curPairT = pairT( ...
    (isempty(curAxonClasses) | ismember(pairT.axonClass, curAxonClasses)) ...
  & (isempty(curTargetClasses) | ismember(pairT.targetClass, curTargetClasses)), :);

[curAxonIds, ~, curPairT.preAggloId] = unique(curPairT.preAggloId);
[curDendIds, ~, curPairT.postAggloId] = unique(curPairT.postAggloId);

curPairT.medLog10AsiArea = cellfun( ...
    @(ids) median(log10(asiT.area(ids))), curPairT.asiIds);

curAreas = [ ...
    curPairT.preAggloId, ...
    curPairT.postAggloId];
curAreas = accumarray( ...
    curAreas, curPairT.medLog10AsiArea, ...
   [numel(curAxonIds), numel(curDendIds)], [], nan);

Util.log('Clustering weight matrix');
[curLink, curDist] = connectEM.Consistency.linkage(curAreas);
Util.log('Done!');
