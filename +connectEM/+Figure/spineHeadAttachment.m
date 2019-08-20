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

shId = 6;
cam = struct;
cam.CameraPosition  = [26.6944, 30.3816, -1.4454];
cam.CameraTarget    = [2.9505, 2.0485, 8.2460];
cam.CameraUpVector  = [-0.6758, 0.6708, 0.3053];
cam.CameraViewAngle = 10.9061;

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

shT.agglo = shAgglos;
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
  & not(shT.preAttached) ...
  & cellfun(@isscalar, shT.dendIds), :);
shT.dendId = cell2mat(shT.dendIds);
shT.dendIds = [];

%% Prototyping
clear cur*;

curSh = shT(shId, :);
curShSegIds = double(curSh.agglo{1});
curDendSegIds = double(conn.dendrites{curSh.dendId});

curNeckSegIds = double(curSh.edges{1}(:, 2));
curNeckSegIds = setdiff(curNeckSegIds, [0; curShSegIds], 'stable');
curNeckSegIds = curNeckSegIds(1:(end - 1));

curDendSegIds = setdiff(curDendSegIds, [curShSegIds; curNeckSegIds]);
curAgglos = [{curShSegIds}; num2cell(curNeckSegIds); {curDendSegIds}];
curLUT = Agglo.buildLUT(maxSegId, curAgglos);

curBox = segPos(curNeckSegIds, :);
curBox = transpose(cat(1, min(curBox, [], 1), max(curBox, [], 1)));
curBox = max(1, curBox - [+256, -256]);

curSeg = loadSegDataGlobal(param.seg, curBox);
curSeg(curSeg ~= 0) = curLUT(curSeg(curSeg ~= 0));
curIsos = arrayfun(@(i) isosurface(curSeg == i, 0.5), 1:numel(curAgglos));

curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');

curColors = prism(numel(curIsos));
curColors = flip(curColors, 1);

for curIdx = 1:numel(curIsos)
    curIso = curIsos(curIdx);
    curIso.vertices = ...
        curIso.vertices ...
     .* voxelSize / 1E3;
    
    curPatch = patch(curIso);
    curPatch.EdgeAlpha = 0;
    curPatch.FaceColor = curColors(curIdx, :);    
end

axis(curAx, 'equal');

curKeys = fieldnames(cam);
for curIdx = 1:numel(curKeys)
    curKey = curKeys{curIdx};
    curAx.(curKey) = cam.(curKey);
end

camlight(curAx);

for curDim = 'XYZ'
    curAxis = sprintf('%sAxis', curDim);
    curAxis = curAx.(curAxis);
    curAxis.Visible = 'off';
end

annotation( ...
    curFig, 'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');
connectEM.Figure.config(curFig, info);
