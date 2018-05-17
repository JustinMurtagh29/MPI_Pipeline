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

%% Calculate spine lengths
spineLengths = ...
    connectEM.Dendrite.calculateSpineLengths( ...
        param, trunks, dendrites, spineHeads);

%% Prepare analysis
maxSegId = Seg.Global.getMaxSegId(param);
spineLUT = Agglo.buildLUT(maxSegId, spineHeads);

synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.type = syn.synapses.type(synT.id);
synT(~synT.isSpine, :) = [];

[~, ~, synT.fold] = unique(synT(:, ...
    {'preAggloId', 'postAggloId'}), 'rows');
fold = accumarray(synT.fold, 1);
synT.fold = fold(synT.fold);

synT.spineId = cellfun( ...
    @(segIds) mode(nonzeros(spineLUT(segIds))), ...
    syn.synapses.postsynId(synT.id));

synT.spineLength = spineLengths(synT.spineId);
synT.spineLength = synT.spineLength / 1E3;
synT = synT(synT.spineLength > 0, :);

synT.preClass = conn.axonMeta.axonClass(synT.preAggloId);
synT.postClass = conn.denMeta.targetClass(synT.postAggloId);

%% Thalamocortical vs. corticocortical
fig = figure();
fig.Color = 'white';

ax = axes(fig);
hold(ax, 'on');

histogram( ...
    ax, synT.spineLength(synT.preClass == 'Thalamocortical'), ...
    'BinEdges', binEdges, 'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);
histogram( ...
    ax, synT.spineLength(synT.preClass == 'Corticocortical'), ...
    'BinEdges', binEdges, 'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);

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

histogram( ...
    ax, synT.spineLength(synT.postClass == 'OtherDendrite'), ...
    'BinEdges', binEdges, 'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);
histogram( ...
    ax, synT.spineLength(synT.postClass == 'ApicalDendrite'), ...
    'BinEdges', binEdges, 'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);
histogram( ...
    ax, synT.spineLength(synT.postClass == 'WholeCell'), ...
    'BinEdges', binEdges, 'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);

leg = legend(ax, ...
    'OtherDendrite', ...
    'ApicalDendrite', ...
    'WholeCell', ...
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
