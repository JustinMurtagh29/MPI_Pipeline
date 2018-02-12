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

%{
synFile = '/home/amotta/Desktop/SynapseAgglos_v3.mat';
connFile = '/home/amotta/Desktop/connectome_axons_18_a_with_den_meta.mat';
%}

synCountFilter = @(n) (n == 2); % couplings with exactly two synapses, only
spineSynFilter = @(s) s; % spine synapses only

%% loading data
% NOTE(amotta): Synapses sizes are contained in the `contactArea` field of 
% `conn.connectomeMeta`. Each cell contains the synapses sizes of the
% correponding entries in `conn.connectome`.
conn = load(connFile);
syn = load(synFile);

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