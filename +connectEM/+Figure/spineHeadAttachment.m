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

shId = 14;

cams = struct;
cams(1).CameraPosition  = [-79.0698, 2.9561, -19.2050];
cams(1).CameraTarget    = [3.1753, 2.9561, 7.5180];
cams(1).CameraUpVector  = [0, 0, 1];
cams(1).CameraViewAngle = 10.9061;

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

%% Rendering
clear cur*;

curSh = shT(shT.id == shId, :);
curShs = shT(shT.dendId == curSh.dendId & shT.autoAttached, :);
curDendSegIds = double(conn.dendrites{curSh.dendId});

curShs.neckSegIds = cellfun( ...
    @(head, neck) setdiff(neck, [0; head], 'stable'), ...
    curShs.agglo, curShs.edges, 'UniformOutput', false);

curSpineSegIds = cell2mat([curShs.agglo; curShs.neckSegIds]);
curDendSegIds = setdiff(curDendSegIds, curSpineSegIds);

curAgglos = cellfun( ...
    @(head, neck) [{head}; num2cell(neck(:))], ...
    curShs.agglo, curShs.neckSegIds, 'UniformOutput', false);

curCmap = jet(120);
curCmap = curCmap((size(curCmap, 1) - 100) / 2 + (1:100), :);

curColors = cellfun( ...
    @(a) feval(@(cols) cols(1:(end - 1), :), curCmap( ...
        ceil(linspace(1, size(curCmap, 1), 1 + numel(a))), :)), ...
    curAgglos, 'UniformOutput', false);

curAgglos = cat(1, curAgglos{:}, {curDendSegIds});
curColors = cat(1, curColors{:}, curCmap(end, :));
curLUT = Agglo.buildLUT(maxSegId, curAgglos);

curBox = segPos(curShs.neckSegIds{curShs.id == shId}, :);
curBox = transpose(cat(1, min(curBox, [], 1), max(curBox, [], 1)));
curBox = max(1, curBox - [+256, -256]);

curSeg = loadSegDataGlobal(param.seg, curBox);
curSeg(curSeg ~= 0) = curLUT(curSeg(curSeg ~= 0));

curIsos = cell(0, 1);
for curIdx = 1:numel(curAgglos)
    curMask = curSeg == curIdx;
    curMask = imclose(curMask, strel('cube', 2));
    curIsos{curIdx} = isosurface(curMask, 0.5);
end

curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');

for curIdx = 1:numel(curIsos)
    curIso = curIsos{curIdx};
    curIso.vertices = reshape(curIso.vertices, [], 3);
    curIso.vertices = curIso.vertices .* voxelSize / 1E3;
    
    curPatch = patch(curIso);
    curPatch.FaceColor = curColors(curIdx, :);
    curPatch.EdgeAlpha = 0;
end

axis(curAx, 'equal');

curKeys = fieldnames(cams);
for curIdx = 1:numel(curKeys)
    curKey = curKeys{curIdx};
    curAx.(curKey) = cams.(curKey);
end

camlight(curAx);

for curDimIdx = 1:3
    curDim = char(double('X') - 1 + curDimIdx);
    curLims = [1, diff(curBox(curDimIdx, :))];
    curLims = curLims * voxelSize(curDimIdx) / 1E3;
    
    curAxis = sprintf('%sAxis', curDim);
    curAxis = curAx.(curAxis);
    curAxis.Visible = 'off';
    curAxis.Limits = curLims;
end

annotation( ...
    curFig, 'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');
connectEM.Figure.config(curFig);
