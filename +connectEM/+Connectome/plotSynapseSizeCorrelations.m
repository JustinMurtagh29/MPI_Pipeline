% TODO(amotta)
% * `contactArea` in `connectomeMeta` ranges up to ~120. That's true even
%   for spine-head filtered synapses. Find out what's going on: Do I
%   misunderstand the data structure or the physical units? Is there a bug?
%   Is synapse agglomeration off? Soma segments classified as spine?
%
% Written by`
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_with_den_meta.mat');

synCountFilter = @(n) (n == 2); % couplings with exactly two synapses, only
spineSynFilter = @(s) s; % spine synapses only

info = Util.runInfo();

%% loading data
% NOTE(amotta): Synapses sizes are contained in the `contactArea` field of 
% `conn.connectomeMeta`. Each cell contains the synapses sizes of the
% correponding entries in `conn.connectome`.
conn = load(connFile);
syn = load(synFile);

% for debugging
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

points = Seg.Global.getSegToPointMap(param);

%% limit synapses
synT = table;
synT.id = cell2mat(conn.connectome.synIdx);
synT.area = cell2mat(conn.connectomeMeta.contactArea);
synT.isSpine = syn.isSpineSyn(synT.id);

synT.preAggloId = repelem( ...
    conn.connectome.edges(:, 1), ...
    cellfun(@numel, conn.connectome.synIdx));

synT.postAggloId = repelem( ...
    conn.connectome.edges(:, 2), ...
    cellfun(@numel, conn.connectome.synIdx));

% limit to spine synapses
synT(~synT.isSpine, :) = [];
synT.isSpine = [];

% remove duplicate entries
[~, uniRows, uniCount] = unique(synT.id);
synT = synT(uniRows, :);

% remove synapses occuring multiple times
% (i.e., between at least two different pairs of neurites)
synT.occurences = accumarray(uniCount, 1);
synT(synT.occurences > 1, :) = [];
synT.occurences = [];

% remove synapses whose size is obviously wrong
synT(synT.area > 1.5, :) = [];

%% look at doubly coupled neurites
[dupNeurites, ~, uniRows] = unique( ...
    synT(:, {'preAggloId', 'postAggloId'}), 'rows');

dupNeurites.areas = accumarray( ...
    uniRows, synT.area, [], ...
    @(a) {reshape(sort(a), 1, [])});

uniCount = accumarray(uniRows, 1);
dupNeurites(uniCount ~= 2, :) = [];

dupNeurites.areas = cell2mat(dupNeurites.areas);

fig = figure();
ax = axes(fig);

hold(ax, 'on');
scatter(ax, ...
    dupNeurites.areas(:, 2), ...
    dupNeurites.areas(:, 1), 12, '+');

xlim([1E-2, 1E1]); xlabel('Axon-spine interface 1 (µm²)');
ylim([1E-2, 1E1]); ylabel('Axon-spine interface 2 (µm²)');

ax.XScale = 'log';
ax.YScale = 'log';

% do the fit
% taken from `matlab/+Analysis/+Script/bartolEtAl2015eLife.m` in `amotta`
xLog = log10(dupNeurites.areas(:, 2));
yLog = log10(dupNeurites.areas(:, 1));

b = [ones(numel(xLog), 1), xLog] \ yLog;
b(1) = 10 ^ b(1);

fitF = @(x) b(1) .* (x .^ b(2));
fitName = sprintf('y = %.2f x^{%.2f}', b(1), b(2));
rawName = sprintf('Raw data (n = %d)', numel(xLog));

fitRange = xlim();
fitRange = linspace(fitRange(1), fitRange(end), 2);
plot(fitRange, fitF(fitRange));
plot(fitRange, fitRange, 'k--');

title({'Same-axon same-dendrite spine synapses'; info.git_repos{1}.hash});
legend(rawName, fitName, 'Location', 'NorthWest');

%% plot distribution of synapse size
% plot distribution
fig = figure();
ax = axes(fig);

histogram(ax, synT.area, linspace(0, 1.5, 51));
xlabel(ax, 'Axon-spine interface (µm²)');
ylabel(ax, 'Spine synapses');
ax.TickDir = 'out';

annotation(...
    fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'String', { ...
        'Spine synapse size distribution';
        info.git_repos{1}.hash});
    
fig.Position(3:4) = [820, 475];

%% write out large synapses for inspection in webKNOSSOS
rng(0);
largeSynT = synT(synT.area > 1.5, :);
largeSynT = largeSynT(randperm(size(largeSynT, 1), 50), :);
numDigits = ceil(log10(1 + size(largeSynT, 1)));

skelNames = arrayfun(@(i, n, a) sprintf( ...
    '%0*d. Synapse %d (%.1f µm²)', numDigits, i, n, a), ...
    1:size(largeSynT, 1), largeSynT.id', largeSynT.area', ...
    'UniformOutput', false);
skelNames = [ ...
    strcat(skelNames, {'. Pre'}); ...
    strcat(skelNames, {'. Post'})];

skelNodes = syn.synapses(largeSynT.id, {'presynId', 'postsynId'});
skelNodes = transpose(table2cell(skelNodes));

skelNodes = cellfun( ...
    @(segIds) points(segIds, :), ...
    skelNodes, 'UniformOutput', false);

skel = Skeleton.fromMST(skelNodes(:), param.raw.voxelSize);
skel.names = skelNames(:);

skel = Skeleton.setParams4Pipeline(skel, param);
skel.write('/home/amotta/Desktop/large-spine-synapses.nml');

%% find doubly coupled neurites
% filter based on spine / non-spine
spineSynMask = spineSynFilter(syn.isSpineSyn);

conn.connectomeMeta.contactArea = cellfun( ...
    @(ids, areas) areas(spineSynMask(ids)), ...
    conn.connectome.synIdx, conn.connectomeMeta.contactArea, ...
    'UniformOutput', false);
conn.connectome.synIdx = cellfun( ...
    @(ids) ids(spineSynMask(ids)), ...
    conn.connectome.synIdx, 'UniformOutput', false);

% filter based on synapse count
conn.connectome.synCount = ...
    cellfun(@numel, conn.connectome.synIdx);
synCountMask = synCountFilter(conn.connectome.synCount);

conn.connectome(~synCountMask, :) = [];
conn.connectomeMeta.contactArea(~synCountMask) = [];

%% plot correlation for pre- & post-coupled synapses
prePostPairSynAreas = conn.connectomeMeta.contactArea;

% sort areas within
prePostPairSynAreas = cellfun( ...
    @(areas) reshape(sort(areas), 1, []), ...
    prePostPairSynAreas, 'UniformOutput', false);

fig = figure();
ax = axes(fig);

data = cell2mat(prePostPairSynAreas);
dataMax = max(data(:));

scatter(ax, data(:, 2), data(:, 1), 18, '+');
xlim(ax, [0, dataMax]);
ylim(ax, [0, dataMax]);