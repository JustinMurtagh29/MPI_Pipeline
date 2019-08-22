% Based on
%   +connectEM/+Number/spineAttachment.m
%   18304841d915e0c7ec2bba008bad5e4c810e29c7
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
autoFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');

% Dendrite to export
dendId = 7;

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = load(connFile);

autoData = load(autoFile);
shEdges = autoData.shEdges;
shAgglos = autoData.shAgglos;
shAutoAttached = autoData.attached;
preSpineFile = autoData.info.param.dendFile;
clear autoData;

preSpineDends = load(preSpineFile);
preSpineDends = preSpineDends.dendrites(preSpineDends.indBigDends);

voxelSize = param.raw.voxelSize;
maxSegId = Seg.Global.getMaxSegId(param);
segPos = Seg.Global.getSegToPointMap(param);

%% Find automatically (but not prematurely) attached spine heads
clear cur*;

preSpineLUT = Agglo.fromSuperAgglo(preSpineDends);
preSpineLUT = Agglo.buildLUT(maxSegId, preSpineLUT);
dendLUT = Agglo.buildLUT(maxSegId, conn.dendrites);

shT = table;
shT.id = reshape(1:numel(shAgglos), [], 1);
shT.agglo = cellfun(@double, shAgglos, 'UniformOutput', false);
shT.edges = shEdges;

shT.autoAttached = shAutoAttached ~= 0;
shT.preAttached = cellfun(@(segIds) ...
    any(preSpineLUT(segIds)), shT.agglo);

shT.dendIds = cellfun( ...
    @(segIds) setdiff(dendLUT(segIds), 0), ...
    shT.agglo, 'UniformOutput', false);
shT.attached = not(cellfun(@isempty, shT.dendIds));

shT = shT( ...
    shT.autoAttached ...
  & cellfun(@isscalar, shT.dendIds), :);
shT.dendId = cell2mat(shT.dendIds);
shT.dendIds = [];

%% Build NML for Heiko
clear cur*;

curShs = shT(shT.dendId == dendId & shT.autoAttached, :);
curDendSegIds = double(conn.dendrites{dendId});

curLin = @(neck) nonzeros(reshape(transpose(neck), [], 1));
curSkipHead = @(head, neck) setdiff(neck, head, 'stable');
curSkipTrunk = @(neck) neck(1:(end - 1));

curShs.neckSegIds = cellfun( ...
    @(head, neck) curSkipHead(head, curSkipTrunk(curLin(neck))), ...
    curShs.agglo, curShs.edges, 'UniformOutput', false);

curSpineSegIds = cell2mat([curShs.agglo; curShs.neckSegIds]);
curDendSegIds = setdiff(curDendSegIds, curSpineSegIds);

curAgglos = cellfun( ...
    @(head, neck) [{head}; num2cell(neck(:))], ...
    curShs.agglo, curShs.neckSegIds, 'UniformOutput', false);

curSkel = skeleton();
curSkel = Skeleton.setParams4Pipeline(curSkel, param);
curSkel = Skeleton.setDescriptionFromRunInfo(curSkel, info);

curMst = @(skel, segIds) ...
    Skeleton.fromMST(segPos(segIds, :), voxelSize, skel);
curNumDigits = ceil(log10(1 + (numel(curAgglos) + 1)));

for curShIdx = 1:numel(curAgglos)
    curShId = curShs.id(curShIdx);
    curShAgglos = curAgglos{curShIdx};
    
    for curNeckIdx = 1:numel(curShAgglos)
        curName = sprintf( ...
            'spine-head-%d_step-%d', ...
            curShId, curNeckIdx - 1);
        curSkel = curMst(curSkel, curShAgglos{curNeckIdx});
        curSkel.names{end} = curName;
    end
end

curSkel = curMst(curSkel, curDendSegIds);
curSkel.names{end} = sprintf('dendrite-%d', dendId);
