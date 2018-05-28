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

calibNml = fullfile( ...
    fileparts(mfilename('fullpath')), 'annotations', ...
    '2012-09-28_ex145_07x2_ROI2017__explorational__amotta__da37a8.nml');

binEdges = linspace(0, 7, 71);

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

trunks = load(trunkFile);
trunks = Agglo.fromSuperAgglo(trunks.dendrites(trunks.indBigDends));
trunks = trunks(Agglo.calculateVolume(param, trunks) > trunkMinSize);

dendrites = load(dendFile);
spineHeads = dendrites.shAgglos;
dendrites = dendrites.dendrites(dendrites.indBigDends);

[conn, syn] = connectEM.Connectome.load(param, connFile);

%% Map synapses onto spine heads (if possible)
maxSegId = Seg.Global.getMaxSegId(param);
spineLUT = Agglo.buildLUT(maxSegId, spineHeads);

syn.synapses.spineId = cellfun( ...
    @(segIds) mode(nonzeros(spineLUT(segIds))), ...
    syn.synapses.postsynId);
    
%% Export random spine heads for ground-truth
segPoints = Seg.Global.getSegToPointMap(param);

rng(0);
randSpineIds = rmmissing(syn.synapses.spineId);
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
        spineHeads(randSpineIds), ...
        'UniformOutput', false), ...
	param.raw.voxelSize, skel);
skel.names = arrayfun( ...
    @(idx, id) sprintf('%0*d. Spine head %d', numDigits, idx, id), ...
	(1:numel(randSpineIds))', randSpineIds, 'UniformOutput', false);
skel.write('/home/amotta/Desktop/random-spine-heads.nml');

%% Calculate spine lengths
spineLengths = ...
    connectEM.Dendrite.calculateSpineLengths( ...
        param, trunks, dendrites, spineHeads);

%% Use manual annotations for calibration
calib = skeleton(calibNml);

calibT = table;
calibT.spineId = regexpi( ...
    calib.names, 'Spine head (\d+)$', 'tokens', 'once');
calibT.spineId = str2double(vertcat(calibT.spineId{:}));
calibT.calibLength = calib.pathLength([], param.raw.voxelSize);
calibT.autoLength = spineLengths(calibT.spineId);

fig = figure();
fig.Color = 'white';

ax = axes(fig);
ax.TickDir = 'out';
axis(ax, 'square');
hold(ax, 'on');

scatter(ax, calibT.calibLength / 1E3, calibT.autoLength / 1E3, 128, '.');

limits = [0, max(ax.XLim(2), ax.YLim(2))];
ax.XLim = limits; ax.YLim = limits;

plot(ax, limits, limits, 'Color', 'black', 'LineStyle', '--');

ticks = union(0, intersect(xticks(ax), yticks(ax)));
xticks(ax, ticks); yticks(ax, ticks);

xlabel(ax, 'True spine length (µm)');
ylabel(ax, 'Calculated spine length (µm)');

title( ...
    ax, {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

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
