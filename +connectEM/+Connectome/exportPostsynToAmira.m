% This script takes the postsynaptic segment equivalence classes contained
% in a connectome, discards the ones with less than `N` synapses, and
% generates isosurfaces.
%
% Note that the postsynaptic segment equivalence classes in the connectome
% are not equal to the dendrite agglomerates. For example, each whole cell
% is split into a soma and a dendrite component.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outputDir = '/tmpscratch/amotta/l4/2018-01-26-postsyn-for-amira';

connFile = fullfile( ...
    rootDir, 'connectomeState', 'connectome_axons_18_a.mat');
minSynCount = 10;

info = Util.runInfo();

%% loading data
fprintf('Loading data... ');
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

conn = load(connFile);
fprintf('done!\n');

%% use WKW for speed-up
param.seg = struct;
param.seg.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
param.seg.backend = 'wkwrap';

%% filter based on synapse count
postSynCount = accumarray( ...
    conn.connectome.edges(:, 2), ...
    cellfun(@numel, conn.connectome.synIdx)', ...
   [numel(conn.dendrites), 1]);

postSynIds = find(postSynCount >= minSynCount);
postSynAgglos = conn.dendrites(postSynIds);

%% generating isosurfaces
Visualization.exportAggloToAmira( ...
    param, postSynAgglos, outputDir, ...
    'reduce', 0.05, 'smoothSizeHalf', 4, 'smoothWidth', 8);

%% writing info
Util.save(fullfile(outputDir, 'info.mat'), info);