% We have isosurfaces for all "large" axons and for dendrites with at least
% `N` synapses. This script generates a MAT file which lists the PLY 
% files of all pre- and postsynaptic entities, respectively.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outputDir = '/home/amotta/Desktop';

preIsoDir = '/tmpscratch/amotta/l4/2018-01-24-axons-18a-isosurfaces';
postIsoDir = '/tmpscratch/amotta/l4/2018-01-26-postsyn-for-amira';

connFile = fullfile( ...
    rootDir, 'connectomeState', 'connectome_axons_18_a.mat');
minSynCount = 10;

info = Util.runInfo();

%% loading data
conn = load(connFile);
[~, connName] = fileparts(connFile);

%% filtering agglomerates
% axons
axonSynCount = conn.axonMeta.synCount;
axonIds = find(axonSynCount >= minSynCount);

% filter postsynaptic targets
postSynCount = accumarray( ...
    conn.connectome.edges(:, 2), ...
    cellfun(@numel, conn.connectome.synIdx)', ...
   [numel(conn.dendrites), 1]);

postIds = find(postSynCount >= minSynCount);

%% build tables
preSyn = table;
preSyn.id = axonIds;
preSyn.plyFile = arrayfun( ...
    @(i) sprintf('iso-%d.ply', i), ...
    axonIds, 'UniformOutput', false);
preSyn.plyFile = fullfile( ...
    preIsoDir, 'ply', preSyn.plyFile);
all(cellfun(@(n) exist(n, 'file'), preSyn.plyFile));

postSyn = table;
postSyn.id = postIds;
postSyn.targetClass = ...
    conn.denMeta.targetClass(postIds);
postSyn.plyFile = arrayfun( ...
    @(i) sprintf('iso-%d.ply', i), ...
    reshape(1:numel(postIds), [], 1), ...
    'UniformOutput', false);
postSyn.plyFile = fullfile( ...
    postIsoDir, 'ply', postSyn.plyFile);

%% write output
outputFile = fullfile(outputDir, sprintf('%s_plys.mat', connName));
Util.save(outputFile, info, preSyn, postSyn);
