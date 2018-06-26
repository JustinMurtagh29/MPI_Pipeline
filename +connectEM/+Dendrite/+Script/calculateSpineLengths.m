% This script calculates the length of spine necks as the shortest path
% from a spine head segment to a trunk segment along the super-agglomerate
% edges.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
trunkFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3.mat');
dendFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto-and-manual.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

% For the origin of this magic constant, see
% https://gitlab.mpcdf.mpg.de/connectomics/amotta/blob/534c026acd534957b84e395d697ac48b3cc6a7ad/matlab/+L4/+Spine/buildDendriteTrunkMask.m
trunkMinSize = 10 ^ 5.5;

calibNml = ...
    connectEM.Dendrite.Data.getFile('spineLengthCalibration.nml');
premCalibNml = ...
    connectEM.Dendrite.Data.getFile('prematureSpineLengthCalibration.nml');

binEdges = linspace(0, 7, 71);

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);

[conn, syn] = connectEM.Connectome.load(param, connFile);

trunks = load(trunkFile);
trunks = Agglo.fromSuperAgglo(trunks.dendrites);
trunks = trunks(Agglo.calculateVolume(param, trunks) > trunkMinSize);
trunkLUT = Agglo.buildLUT(maxSegId, trunks);

dendrites = load(dendFile);
shAgglos = dendrites.shAgglos;
dendrites = dendrites.dendrites;

% Only keep dendrites that also contain a trunk
dendMask = Agglo.fromSuperAgglo(dendrites, true);
dendMask = cellfun(@(ids) any(trunkLUT(ids)), dendMask);
dendrites = dendrites(dendMask);

% Sanity check
assert(numel(trunks) == numel(dendrites));

%% Build spine head table
shT = table;
shT.agglo = shAgglos;

dendLUT = Agglo.fromSuperAgglo(dendrites, true);
dendLUT = Agglo.buildLUT(maxSegId, dendLUT);

shT.attached = cellfun(@(ids) any(dendLUT(ids)), shT.agglo);
shT.premAttached = cellfun(@(ids) any(trunkLUT(ids)), shT.agglo);

%% Map synapses onto spine heads (if possible)
spineLUT = Agglo.buildLUT(maxSegId, shT.agglo);
trunkLUT = Agglo.buildLUT(maxSegId, trunks);

syn.synapses.spineId = cellfun( ...
    @(segIds) mode(nonzeros(spineLUT(segIds))), ...
    syn.synapses.postsynId);
syn.synapses.trunkId = cellfun( ...
    @(segIds) mode(nonzeros(trunkLUT(segIds))), ...
    syn.synapses.postsynId);
    
%% Export random spine heads for ground-truth
rng(0);

%{
% Random spine heads
randSpineIds = rmmissing(syn.synapses.spineId);
%}

% Random prematurely attached spines
% I.e., spines attached before running the spine attachment routines
randSpineIds = isnan(syn.synapses.spineId) | isnan(syn.synapses.trunkId);
randSpineIds = syn.synapses.spineId(~randSpineIds);

% Select random subset
randSpineIds = randSpineIds(randperm(numel(randSpineIds)));
randSpineIds = randSpineIds(1:25);

numDigits = ceil(log10(1 + numel(randSpineIds)));

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));
skel = Skeleton.fromMST( ...
    cellfun( ...
        @(segIds) segPoints(segIds, :), ...
        shT.agglo(randSpineIds), ...
        'UniformOutput', false), ...
	param.raw.voxelSize, skel);
skel.names = arrayfun( ...
    @(idx, id) sprintf('%0*d. Spine head %d', numDigits, idx, id), ...
	(1:numel(randSpineIds))', randSpineIds, 'UniformOutput', false);
skel.write('/home/amotta/Desktop/random-spine-heads.nml');

%% Calculate spine lengths
shT.length = ...
    connectEM.Dendrite.calculateSpineLengths( ...
        param, trunks, dendrites, shT.agglo);

%% Calibration of automatically or manually attached spine heads
% Note that this excludes "prematurely" and non-attached spines
calibT = loadCalibrationNml(param, calibNml);
calibT.autoLength = shT.length(calibT.spineId);

calibT(isnan(calibT.autoLength), :) = [];
calibT(~calibT.autoLength, :) = [];

corrCoeff = ...
    sum(calibT.calibLength) ...
  / sum(calibT.autoLength);

premCalibT = loadCalibrationNml(param, premCalibNml);
premCorr = mean(premCalibT.calibLength);

fig = figure();
fig.Color = 'white';

ax = axes(fig);
ax.TickDir = 'out';
axis(ax, 'square');
hold(ax, 'on');

scatter(ax, ...
    calibT.autoLength / 1E3, ...
    calibT.calibLength / 1E3, ...
    128, '.');

limits = [0, max(ax.XLim(2), ax.YLim(2))];
ax.XLim = limits; ax.YLim = limits;

plot(ax, limits, limits, 'Color', 'black', 'LineStyle', '--');
plot(ax, limits, corrCoeff .* limits);

ticks = union(0, intersect(xticks(ax), yticks(ax)));
xticks(ax, ticks); yticks(ax, ticks);

xlabel(ax, 'Shortest path-based spine length (µm)');
ylabel(ax, 'Tracing-based spine length (µm)');

leg = legend(ax, { ...
    sprintf('%d calibration points', height(calibT)), ...
    sprintf('Correction coefficient: %.3g', corrCoeff)}, ...
    'Location', 'NorthWest');
leg.Box = 'off';

title( ...
    ax, {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Estimate total spine length
corrLength = corrCoeff * sum(shT.length(shT.attached & ~shT.premAttached));
premLength = premCorr * sum(shT.attached & shT.premAttached);

totalLenght = (corrLength + premLength) / 1E9;
fprintf('Total spine neck length: %g m\n', totalLenght);

%% Prepare analysis
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.spineId = syn.synapses.spineId(synT.id);
synT.type = syn.synapses.type(synT.id);
synT(~synT.isSpine, :) = [];

[~, ~, synT.fold] = unique(synT(:, ...
    {'preAggloId', 'postAggloId'}), 'rows');
fold = accumarray(synT.fold, 1);
synT.fold = fold(synT.fold);

synT.spineLength = spineLengths(synT.spineId);
synT.spineLength = synT.spineLength / 1E3;

synT.preClass = conn.axonMeta.axonClass(synT.preAggloId);
synT.postClass = conn.denMeta.targetClass(synT.postAggloId);

%% Check if prematurely attached spines are uniformly distributed
premFracPreClasses = accumarray( ...
    double(synT.preClass), synT.spineLength == 0, [], @mean);
premFracPostClasses = accumarray( ...
    double(synT.postClass), synT.spineLength == 0, [], @mean);

%% Clean-up
% Get rid of spines which were already attached prior to running the spine
% attachment routines. For these we cannot calculate a spine length.
% synT = synT(synT.spineLength > 0, :);
synT(isnan(synT.spineLength), :) = [];

%% Thalamocortical vs. corticocortical
curTcData = synT.spineLength(synT.preClass == 'Thalamocortical');
curCcData = synT.spineLength(synT.preClass == 'Corticocortical');

fig = figure();
fig.Color = 'white';

ax = axes(fig);
hold(ax, 'on');

histogram( ...
    ax, curTcData, ...
    'BinEdges', binEdges, 'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);
histogram( ...
    ax, curCcData, ...
    'BinEdges', binEdges, 'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);
plot( ...
    ax, repelem(mean(curTcData), 1, 2), ...
    ax.YLim, 'Color', ax.ColorOrder(1, :), 'LineStyle', '--');
plot( ...
    ax, repelem(mean(curCcData), 1, 2), ...
    ax.YLim, 'Color', ax.ColorOrder(2, :), 'LineStyle', '--');

leg = legend(ax, ...
    'Thalamocortical', ...
    'Corticocortical', ...
    'Location', 'NorthEast');
leg.Box = 'off';

ax.TickDir = 'out';
ax.XLim = binEdges([1, end]);
xlabel(ax, 'Spine length (µm)');
ylabel(ax, 'Probability');

title( ...
    ax, {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Additional inhibitory input
spineIds = unique(synT.spineId);
spineHasInh = ismember(spineIds, ...
    synT.spineId(synT.type == 'SecondarySpine'));

fig = figure();
fig.Color = 'white';

ax = axes(fig);
hold(ax, 'on');

histogram( ...
    ax, spineLengths(spineIds(spineHasInh)) / 1E3, ...
    'BinEdges', binEdges, 'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);
histogram( ...
    ax, spineLengths(spineIds(~spineHasInh)) / 1E3, ...
    'BinEdges', binEdges, 'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);

leg = legend(ax, ...
    'With inhibition', ...
    'Without inhibition', ...
    'Location', 'NorthEast');
leg.Box = 'off';

ax.TickDir = 'out';
ax.XLim = binEdges([1, end]);
xlabel(ax, 'Spine length (µm)');
ylabel(ax, 'Probability');

title( ...
    ax, {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Slightly off-topic analysis
% Are TC spine synapses more likely to be co-innervated by an IN axon?
curSynT = sortrows(synT, {'spineId', 'type'});
[~, uniRows] = unique(curSynT.spineId);
curSynT = curSynT(uniRows, :);

curSynT.hasInh = ismember(curSynT.spineId, ...
    synT.spineId(synT.type == 'SecondarySpine'));

curResT = table;
curResT.synType = unique(curSynT.preClass);
curResT.inhFrac = accumarray( ...
    double(curSynT.preClass), curSynT.hasInh, ...
    size(curResT.synType), @mean);

disp(curResT);

%% Dendrite type
fig = figure();
fig.Color = 'white';

ax = axes(fig);
hold(ax, 'on');

curOdData = synT.spineLength(synT.postClass == 'OtherDendrite');
curAdData = synT.spineLength(synT.postClass == 'ApicalDendrite');

histogram( ...
    ax, curOdData, ...
    'BinEdges', binEdges, 'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);
histogram( ...
    ax, curAdData, ...
    'BinEdges', binEdges, 'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);
plot( ...
    ax, repelem(mean(curOdData), 1, 2), ...
    ax.YLim, 'Color', ax.ColorOrder(1, :), 'LineStyle', '--');
plot( ...
    ax, repelem(mean(curAdData), 1, 2), ...
    ax.YLim, 'Color', ax.ColorOrder(2, :), 'LineStyle', '--');
leg = legend(ax, ...
    'OtherDendrite', ...
    'ApicalDendrite', ...
    'Location', 'NorthEast');
leg.Box = 'off';

ax.TickDir = 'out';
ax.XLim = binEdges([1, end]);
xlabel(ax, 'Spine length (µm)');
ylabel(ax, 'Probability');

title( ...
    ax, {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Variability across soma-based reconstructions
curBinEdges = linspace(0, 3, 31);
curSynT = synT(synT.postClass == 'WholeCell', :);

curCellT = table;
[curCellT.postAggloId, ~, curCellSpineLengths] = unique(curSynT.postAggloId);
curCellT.id = conn.denMeta.cellId(curCellT.postAggloId);
curCellT.isInterneuron = conn.denMeta.isInterneuron(curCellT.postAggloId);

curCellT.spineLengths = accumarray( ...
    curCellSpineLengths, curSynT.spineLength, [], @(vals) {vals});
curCellT.meanSpineLength = cellfun(@mean, curCellT.spineLengths);

curCellT(cellfun(@numel, curCellT.spineLengths) < 50, :) = [];

fig = figure;
fig.Color = 'white';

ax = axes(fig);
hold(ax, 'on');

histogram( ...
    ax, curCellT.meanSpineLength(~curCellT.isInterneuron), ...
    'BinEdges', curBinEdges, 'DisplayStyle', 'stairs', 'LineWidth', 2);
histogram( ...
    ax, curCellT.meanSpineLength(curCellT.isInterneuron), ...
    'BinEdges', curBinEdges, 'DisplayStyle', 'stairs', 'LineWidth', 2);

ax.TickDir = 'out';
ax.XLim = curBinEdges([1, end]);
xlabel(ax, 'Mean spine length (µm)');
ylabel(ax, 'Neurons');

leg = legend( ...
    ax, 'Excitatory neurons', 'Interneurons', ...
    'Location', 'NorthWest');
leg.Box = 'off';

title( ...
    ax, {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Synapse size
fig = figure();
fig.Color = 'white';

ax = axes(fig);
hold(ax, 'on');

fitRes = fit(synT.spineLength, synT.area, 'poly1');
scatter(ax, synT.spineLength, synT.area, '.');
plot(ax, xlim(ax), fitRes(xlim(ax)), 'Color', 'black');

ax.TickDir = 'out';
xlabel(ax, 'Spine length (µm)');
ylabel(ax, 'Axon-spine interface area (µm²)');

leg = legend(ax, ...
    sprintf('Data points (n = %d)', height(synT)), ...
    sprintf('Linear fit (y = %f + %fx)', fitRes.p2, fitRes.p1), ...
    'Location', 'NorthEast');
leg.Box = 'off';

title( ...
    ax, {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Number of synapses between axon-dendrite pair
plotFolds = 1:4;

fig = figure();
fig.Color = 'white';

ax = axes(fig);
hold(ax, 'on');

for curFold = plotFolds
    histogram( ...
        ax, synT.spineLength(synT.fold == curFold), ...
        'BinEdges', binEdges, 'Normalization', 'probability', ...
        'DisplayStyle', 'stairs', 'LineWidth', 2);
end

% Plot mean values
for curFoldIdx = 1:numel(plotFolds)
    curFold = plotFolds(curFoldIdx);
    curMean = mean(synT.spineLength(synT.fold == curFold));
    
    plot( ...
        ax, [curMean, curMean], ax.YLim, ...
        'Color', ax.ColorOrder(curFoldIdx, :), ...
        'LineStyle', '--');
end

ax.TickDir = 'out';
ax.XLim = binEdges([1, end]);
xlabel(ax, 'Spine length (µm)');
ylabel(ax, 'Axon-spine interface area (µm²)');

legends = arrayfun( ...
    @(n) sprintf('%d joint spine synapse(s)', n), ...
    plotFolds, 'UniformOutput', false);
leg = legend(ax, legends, 'Location', 'NorthEast');
leg.Box = 'off';

title( ...
    ax, {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Check for correlation between synapse size and spine neck length
% among joint spine synapses.
[~, ~, synT.jointId] = unique(synT(:, ...
    {'preAggloId', 'postAggloId'}), 'rows');

jointFold = accumarray(synT.jointId, 1);
jointAreas = accumarray(synT.jointId, synT.area, [], @(vals) {vals});
jointLengths = accumarray(synT.jointId, synT.spineLength, [], @(vals) {vals});

corrT = table;
corrT.fold = transpose(2:4);
corrT.rankCorr = nan(size(corrT.fold));

for curIdx = 1:height(corrT)
    curFold = corrT.fold(curIdx);
    curJointAreas = jointAreas(jointFold == curFold);
    curJointLengths = jointLengths(jointFold == curFold);

    curJointCorr = cellfun( ...
        @(a, l) corr(a(:), l(:), 'Type', 'Kendall'), ...
        curJointAreas, curJointLengths);
    corrT.rankCorr(curIdx) = mean(curJointCorr, 'omitnan');
end

%% Utilities
function calibT = loadCalibrationNml(param, nmlFile)
    calib = skeleton(nmlFile);

    calibT = table;
    calibT.spineId = regexpi( ...
        calib.names, 'Spine head (\d+)$', 'tokens', 'once');
    calibT.spineId = str2double(vertcat(calibT.spineId{:}));
    calibT.calibLength = calib.pathLength([], param.raw.voxelSize);
end
