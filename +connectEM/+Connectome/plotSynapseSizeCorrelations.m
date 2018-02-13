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

%% plot distribution of synapse size
synT = table;
synT.id = cell2mat(conn.connectome.synIdx);
synT.area = cell2mat(conn.connectomeMeta.contactArea);

% limit to spine synapses
synT(~syn.isSpineSyn(synT.id), :) = [];

% remove duplicate entries
[~, uniRows] = unique(synT.id);
synT = synT(uniRows, :);

fprintf('Spine synapses\n');
fprintf('  %d spine synapses in connectome\n', size(synT, 1));
fprintf('  %d spine synapses with ASI > 1.5 µm²\n', sum(synT.area > 1.5));
fprintf('  Largest ASI: %.1f µm\n', max(synT.area));

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
        sprintf('%s, figure %d', mfilename, fig.Number);
        info.git_repos{1}.hash});
    
fig.Position(3:4) = [820, 475];

%% write out large synapses for inspection in webKNOSSOS
rng(0);
largeSynT = synT(synT.area > 1.5, :);
largeSynT = largeSynT(randperm(size(largeSynT, 1), 50), :);
numDigits = ceil(log10(1 + size(largeSynT, 1)));

skelNames = arrayfun(@(i, n, a) sprintf( ...
    '%0*d. Synapse %d (%.1f µm²)', numDigits, i, n, a), ...
    1:size(synIds, 1), largeSynT.id', largeSynT.area', ...
    'UniformOutput', false);
skelNames = [ ...
    strcat(skelNames, {'. Pre'}); ...
    strcat(skelNames, {'. Post'})];

skelNodes = syn.synapses(synIds, {'presynId', 'postsynId'});
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