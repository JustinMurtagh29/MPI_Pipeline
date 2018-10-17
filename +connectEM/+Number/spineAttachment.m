% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
autoFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = load(connFile);

autoData = load(autoFile);
shAgglos = autoData.shAgglos;
shAutoAttached = autoData.attached;
preSpineFile = autoData.info.param.dendFile;
clear autoData;

preSpineDends = load(preSpineFile);
preSpineDends = preSpineDends.dendrites(preSpineDends.indBigDends);

box = param.bbox;
voxelSize = param.raw.voxelSize;
maxSegId = Seg.Global.getMaxSegId(param);

segWeights =  Seg.Global.getSegToSizeMap(param);
segPoints = Seg.Global.getSegToCentroidMap(param);

%% Process spine heads
preSpineLUT = Agglo.fromSuperAgglo(preSpineDends);
preSpineLUT = Agglo.buildLUT(maxSegId, preSpineLUT);
dendLUT = Agglo.buildLUT(maxSegId, conn.dendrites);

shT = table;
shT.agglo = shAgglos;
shT.autoAttached = shAutoAttached ~= 0;
shT.preAttached = cellfun(@(segIds) any(preSpineLUT(segIds)), shT.agglo);
shT.attached = cellfun(@(segIds) any(dendLUT(segIds)), shT.agglo);

wmean = @(w, v) sum(v .* (w / sum(w)), 1);
shT.pos = cell2mat(cellfun(@(segIds) ...
    wmean(segWeights(segIds), segPoints(segIds, :)), ...
    shT.agglo, 'UniformOutput', false));

shT.borderDist = min( ...
    shT.pos - transpose(box(:, 1)), ...
    transpose(box(:, 2)) - shT.pos);
shT.borderDist = shT.borderDist .* voxelSize;
shT.borderDist = min(shT.borderDist, [], 2);

%% Prepare plot
shT = sortrows(shT, 'borderDist', 'descend');

countVsDist = reshape(1:height(shT), [], 1);
attachedVsDist = cumsum(shT.attached);
borderDist = shT.borderDist;

fracAttachedVsDist = attachedVsDist ./ countVsDist;

%% Plot
curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [355, 345];

curAx = axes(curFig);
plot(curAx, borderDist / 1E3, fracAttachedVsDist, 'LineWidth', 2);
axis(curAx, 'square');

curAx.Box = 'off';
curAx.TickDir = 'out';
curAx.XLim = [0, borderDist(1) / 1E3];
curAx.YLim = [0, 1];

xlabel(curAx, 'Minimum distance from dataset border (µm)');
ylabel(curAx, 'Fraction of spine heads attached');

title(curAx, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Quantitative evaluation
numSpineHeads = height(shT) %#ok
numManuallyAttachedSpineHeads = sum(shT.attached & ~shT.autoAttached) %#ok

fractionOfSpineHeadsPrematurelyAttached = mean(shT.preAttached) %#ok
fractionOfSpineHeadsPrematurelyAttachedAfterTenUm = ...
    mean(shT.preAttached(shT.borderDist > 10E3)) %#ok

fractionOfSpineHeadsAutoAttached = mean(shT.autoAttached) %#ok
fractionOfSpineHeadsAutoAttachedAfterTenUm = ...
    mean(shT.autoAttached(shT.borderDist > 10E3)) %#ok

fracionOfSpineHeadsAttached = mean(shT.attached) %#ok
fractionOfSpineHeadsAttachedAfterTenUm = ...
    mean(shT.attached(shT.borderDist > 10E3)) %#ok
