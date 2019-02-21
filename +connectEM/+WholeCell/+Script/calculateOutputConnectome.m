% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');

% NOTE(amotta): This file is identical to 20180726T190355_results.mat with
% the exception of the additional `outputMap.axonData.segIds` field.
outputMapFile = '/tmpscratch/amotta/l4/2018-07-26-tracing-based-output-maps/20190117T143833_results.mat';

runId = '20190221T112510';

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);

shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

graph = Graph.load(rootDir);
graph = graph(graph.borderIdx ~= 0, :);

% Loading tracings and synapses of L4 cells
outputMap = load(outputMapFile);
connFile = outputMap.info.param.connFile;
axonData = outputMap.axonData;

% Load connectome with synapses and TC axons
[conn, syn] = connectEM.Connectome.load(param, connFile);

%%
% NOTE(amotta): Because the axons used for this analysis are different
% from the ones in the dense reconstruction, we have to build a fake
% connectome.
clear cur*;

axonAgglos = reshape({axonData.segIds}, [], 1);
axonAgglos = Agglo.removeSegmentOverlap(axonAgglos);

axonLUT = Agglo.buildLUT(maxSegId, axonAgglos);
dendLUT = Agglo.buildLUT(maxSegId, conn.dendrites);

curSyn = syn.synapses;
curSyn.id = reshape(1:height(curSyn), [], 1);

curSyn.preAggloId = cellfun( ...
    @(ids) setdiff(axonLUT(ids), 0), ...
    curSyn.presynId, 'UniformOutput', false);
curSyn = curSyn(cellfun( ...
    @isscalar, curSyn.preAggloId), :);

curSyn.postAggloId = cellfun( ...
    @(ids) setdiff(dendLUT(ids), 0), ...
    curSyn.postsynId, 'UniformOutput', false);
curSyn = curSyn(cellfun( ...
    @isscalar, curSyn.postAggloId), :);

curSyn = [ ...
    cell2mat(curSyn.preAggloId), ...
    cell2mat(curSyn.postAggloId), ...
    curSyn.id];

curConnectome = table;
[curConnectome.edges, ~, curIds] = ...
    unique(curSyn(:, 1:2), 'rows');
curConnectome.synIdx = accumarray( ...
    curIds, curSyn(:, end), [], @(ids) {ids});

l4 = struct;
l4.conn = struct;
l4.conn.axons = axonAgglos;
l4.conn.denMeta = conn.denMeta;
l4.conn.dendrites = conn.dendrites;
l4.conn.connectome = curConnectome;

l4.synT = ...
    connectEM.Connectome.buildSynapseTable(l4.conn, syn);
l4.asiT = ...
    connectEM.Connectome.buildAxonSpineInterfaces( ...
        param, graph, shAgglos, l4.conn, syn, 'addBorderIdsVar', true);
    
%% Calculate position
clear cur*;

% Calculate position
curBorderMeta = fullfile(rootDir, 'globalBorder.mat');
curBorderMeta = load(curBorderMeta, 'borderSize', 'borderCoM');
curBorderMeta = structfun(@double, curBorderMeta, 'UniformOutput', false);

% See also connectEM.Consistency.Script.buildAxonSpineInterfaces
curWmean = @(w, v) sum((w / sum(w, 1)) .* v, 1);

l4.asiT.pos = cellfun( ...
    @(ids) curWmean( ...
        curBorderMeta.borderSize(ids), ...
        curBorderMeta.borderCoM(ids, :)), ...
    l4.asiT.borderIds, 'UniformOutput', false);

l4.asiT.pos = round(cell2mat(l4.asiT.pos));
l4.asiT.pos(cellfun(@isempty, l4.asiT.borderIds), :) = nan;

%% Calculate ASI areas
clear cur*;

curIds = find(not(any(isnan(l4.asiT.pos), 2)));
curAxonAgglos = {axonData.segIds};

curAreas = nan(numel(curIds), 1);
parfor curIdx = 1:numel(curIds)
    curId = curIds(curIdx);

    curAreas(curIdx) = ...
        connectEM.Consistency.buildAxonSpineInterfaceAreas( ...
            param, curAxonAgglos, shAgglos, l4.asiT(curId, :)); %#ok
end

l4.asiT.area = nan(height(l4.asiT), 1);
l4.asiT.area(curIds) = curAreas;

%% Save output
clear cur*;

[curDir, curName] = fileparts(outputMapFile);
outFile = sprintf('%s__%s_connectome.mat', curName, runId);
outFile = fullfile(curDir, outFile);

out = l4;
out.info = info;
Util.saveStruct(outFile, out);
Util.protect(outFile);
