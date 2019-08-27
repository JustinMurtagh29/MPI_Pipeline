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

% HACKHACKHACK(amotta): The above connectome was generated based on axons
% 18a. There was, however, a bug in the routine for merging super-agglos
% on overlap which introduced bogus edges. This bug was later fixed,
% resulting in axons 18b. The two super-agglo states 18a and 18b are
% (almost) identical in terms of segment equivalence classes. We need
% proper edges here. So, let's use axons 18b.
axonFile = fullfile(rootDir, 'aggloState', 'axons_18_b.mat');

minSynCount = 10;

info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = load(connFile);
[~, connName] = fileparts(connFile);

axons = load(axonFile);
axons = axons.axons(axons.indBigAxons);

% HACKHACKHACK(amotta): Early versions of the connectome generator did not
% keep track of the mapping from super-agglomerates to segment equivalence
% classes as they are stored in the connectome. Now we have to pull some
% dirty hacks to work around that.
maxSegId = Seg.Global.getMaxSegId(param);
curLUT = Agglo.fromSuperAgglo(axons, true);
curLUT = Agglo.buildLUT(maxSegId, curLUT);

curParentIds = cellfun( ...
    @(segIds) setdiff(curLUT(segIds), 0), ...
    conn.axons, 'UniformOutput', false);

% HACKHACKHACK(amotta): It's possible that a super-agglomerate consisting
% of nothing more than a long flight path reaches the 5 µm size threshold.
% If this super-agglomerate only gets assigned segments during pick-up in
% the connectome generation, there is no way to reconstruct the parent ID
% here. Let's just set the parent ID to zero here and deal with it below.
curParentIds(cellfun(@isempty, curParentIds)) = {0};
curParentIds = cell2mat(curParentIds);

conn.axonMeta.parentId = curParentIds(:);

%% filtering agglomerates
clear cur*;

% axons
axonSynCount = conn.axonMeta.synCount;
axonIds = find(axonSynCount >= minSynCount);

% filter postsynaptic targets
postSynCount = accumarray( ...
    conn.connectome.edges(:, 2), ...
    cellfun(@numel, conn.connectome.synIdx)', ...
   [numel(conn.dendrites), 1]);

postIds = find(postSynCount >= minSynCount);

%% build skeletons
clear cur*;

curEmtpySkel = skeleton();
curEmptySkel = Skeleton.setDescriptionFromRunInfo(curEmtpySkel, info);
curEmptySkel = Skeleton.setParams4Pipeline(curEmptySkel, param);

curOutDir = fullfile(preIsoDir, 'nml');
if ~exist(curOutDir, 'dir'); mkdir(curOutDir); end

for curAxonId = 1:height(conn.axonMeta)
    curParentId = conn.axonMeta.parentId(curAxonId);
    if ~curParentId; continue; end
    
    curSkel = axons(curParentId);
    curSkel = Superagglos.buildAggloAndFlightSkels(curSkel, curEmptySkel);
    
    curDesc = sprintf('Axon %d', curAxonId);
    curSkel = curSkel.setDescription(curDesc, 'append', true);
    
    curOutFile = sprintf('axon-%d.nml', curAxonId);
    curOutFile = fullfile(curOutDir, curOutFile);
    curSkel.write(curOutFile);
end

Util.protect(curOutDir, true);

%% build tables
preSyn = table;
preSyn.id = axonIfds;
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
