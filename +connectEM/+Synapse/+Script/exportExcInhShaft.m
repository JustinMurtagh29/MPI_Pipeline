% This script export random example of shaft synapses made by
% * excitatory (mainly spine targeting) axons, and
% * inhibitory, soma-targeting axons.
%
% It also generates NML files of non-soma targeting, shaft preferring
% axons. This might be an excitatory axon subpopulation from L6.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
outputDir = '/home/amotta/Desktop';

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

points = Seg.Global.getSegToPointMap(param);

[conn, syn, axonClasses] = connectEM.Connectome.load(param, connFile);

axons = load(conn.info.param.axonFile);
axons = axons.axons;

%% Select synapses
% Find non-soma shaft synapses
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.isSoma = conn.denMeta.targetClass(synT.postAggloId) == 'Somata';
synT(synT.isSpine | synT.isSoma, :) = [];

% Excitatory shaft synapses
excAxonIds = axonClasses(1).axonIds;

% Shaft synapses of soma-targeting inhibitory neurons
inhAxonIds = conn.denMeta.id(conn.denMeta.targetClass == 'Somata');
inhAxonIds = ismember(conn.connectome.edges(:, 2), inhAxonIds);
inhAxonIds = unique(conn.connectome.edges(inhAxonIds, 1));
inhAxonIds = intersect(axonClasses(2).axonIds, inhAxonIds);

% Shaft synapses of non-soma-targeting, shaft preferring axon
excShaftAxonIds = setdiff(axonClasses(2).axonIds, inhAxonIds);

exportSyns = struct;
exportSyns(1).tag = 'exc-shaft-synapses';
exportSyns(1).syns = synT(ismember(synT.preAggloId, excAxonIds), :);
exportSyns(2).tag = 'inh-shaft-synapses';
exportSyns(2).syns = synT(ismember(synT.preAggloId, inhAxonIds), :);

%% Export random examples
skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

for curSyn = reshape(exportSyns, 1, [])
    curTag = curSyn.tag;
    curSyns = curSyn.syns;
    
    rng(0);
    curRandIds = randperm(size(curSyns, 1));
    curSyns = curSyns(curRandIds(1:50), :);
    
    curSynAgglos = cellfun( ...
        @vertcat, ...
        syn.synapses.presynId(curSyns.id), ...
        syn.synapses.postsynId(curSyns.id), ...
        'UniformOutput', false);
    curSynAgglos = cellfun( ...
        @(segIds) points(segIds, :), ...
        curSynAgglos, 'UniformOutput', false);
    
    curSkel = Skeleton.fromMST( ...
        curSynAgglos, param.raw.voxelSize, skel);
    
    curNumDigits = ceil(log10(1 + size(curSyns, 1)));
    curSkel.names = arrayfun( ...
        @(idx, id) sprintf('%0*d. Synapse %d', curNumDigits, idx, id), ...
        reshape(1:size(curSyns, 1), [], 1), curSyns.id, ...
        'UniformOutput', false);
    
    curFile = sprintf('%s.nml', curTag);
    curSkel.write(fullfile(outputDir, curFile));
end

%% Export samples of non-soma innervating, shaft preferring axons
rng(0);
curRandIds = randperm(numel(excShaftAxonIds));
curRandIds = excShaftAxonIds(curRandIds);
curRandIds = curRandIds(1:20);

numDigits = ceil(log10(1 + numel(curRandIds)));

for curIdx = 1:numel(curRadIds)
    curId = curRandIds(curIdx);
    curAgglo = axons(conn.axonMeta.parentId(curId));
    
    curSynIds = synT.preAggloId == curId;
    curSynIds = synT.id(curSynIds);
    
    curSynAgglos = cellfun( ...
        @vertcat, ...
        syn.synapses.presynId(curSynIds), ...
        syn.synapses.postsynId(curSynIds), ...
        'UniformOutput', false);
    curSynAgglos = cellfun( ...
        @(segIds) points(segIds, :), ...
        curSynAgglos, 'UniformOutput', false);
    
    curSkel = Superagglos.toSkel(curAgglo, skel);
    curSkel.names{1} = sprintf('Axon %d', curId);
    
    curSkel = Skeleton.fromMST( ...
        curSynAgglos, param.raw.voxelSize, curSkel);
    curSkel.names(2:end) = arrayfun( ...
        @(id) sprintf('Synapse %d', id), ...
        curSynIds, 'UniformOutput', false);
    curSkel.colors(2:end) = {[0, 0, 1, 1]};
    
    curFile = sprintf( ...
        '%0*d_non-soma-shaft-axons_axon-%d.nml', ...
        numDigits, curIdx, curId);
    curSkel.write(fullfile(outputDir, curFile));
end
