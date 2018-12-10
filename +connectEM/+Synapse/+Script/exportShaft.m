% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
outputDir = '/home/amotta/Desktop';

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

points = Seg.Global.getSegToPointMap(param);

[conn, syn] = connectEM.Connectome.load(param, connFile);

%% Select synapses
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT(synT.isSpine, :) = [];

targetClasses = unique(conn.denMeta.targetClass);
targetSynT = arrayfun(@(t) ...
    synT(ismember(synT.postAggloId, ...
        find(conn.denMeta.targetClass == t)), :), ...
    targetClasses, 'UniformOutput', false);

%% Export random examples
exportSyns = lower(arrayfun(@char, targetClasses, 'UniformOutput', false));
exportSyns = struct('tag', exportSyns(:), 'syns', targetSynT(:));

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
    
    curFile = sprintf('shaft-synapses_onto-%s.nml', curTag);
    curSkel.write(fullfile(outputDir, curFile));
end
