% This script export random example of shaft synapses made by
% * excitatory (mainly spine targeting) axons
% * inhibitory, soma-targeting axons
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
