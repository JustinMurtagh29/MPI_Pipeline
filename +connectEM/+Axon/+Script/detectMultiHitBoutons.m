% Detection of multi-hit / -synaptic boutons.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

synPairCount = 100;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);

[conn, syn] = connectEM.Connectome.load(param, connFile);

axons = load(conn.info.param.axonFile, 'axons');
axons = axons.axons;

[connDir, connName] = fileparts(connFile);
interSynFile = sprintf('%s_intersynapse.mat', connName);
interSynFile = fullfile(connDir, interSynFile);
interSyn = load(interSynFile);

%% Training data
% This section select N random synapses, and exports it together with one
% of its direct neighbours. I'll then decide whether or not a pair of
% synapses was originating from the same bouton.
interSyn.synCounts = cellfun(@numel, interSyn.synIds);
randAxonIdx = find(interSyn.synCounts > 1);

rng(0);
randAxonIdx = datasample( ...
    randAxonIdx, synPairCount, ...
    'weights', interSyn.synCounts(randAxonIdx));
randAxonIds = interSyn.axonIds(randAxonIdx);

randSynIds = nan(synPairCount, 2);
for curIdx = 1:synPairCount
    curAxonIdx = randAxonIdx(curIdx);
    
    curSynIds = interSyn.synIds{curAxonIdx};
    curSynToSynDists = interSyn.synToSynDists{curAxonIdx};
    curSynToSynDists(isinf(curSynToSynDists)) = 0;
    
    % Linearize axon
   [~, curRefIdx] = max(curSynToSynDists(:));
   [~, curRefIdx] = ind2sub(size(curSynToSynDists), curRefIdx);
    curRefDists = curSynToSynDists(:, curRefIdx);
    
    % Pick random synapse
    curSynIdx = randi(numel(curSynIds));
    curSynDists = curSynToSynDists(:, curSynIdx);
    
    % Pick random neighbour
    curLeftIdx = setdiff(find( ...
        curRefDists <= curRefDists(curSynIdx)), curSynIdx);
    curRightIdx = setdiff(find( ...
        curRefDists >= curRefDists(curSynIdx)), curSynIdx);
    
   [~, curTempIdx] = min(curSynDists(curLeftIdx));
    curLeftIdx = curLeftIdx(curTempIdx);
    
   [~, curTempIdx] = min(curSynDists(curRightIdx));
    curRightIdx = curRightIdx(curTempIdx);
    
    curNeighIdx = [curLeftIdx, curRightIdx];
    curNeighIdx = curNeighIdx(randi(numel(curNeighIdx)));
    
    randSynIds(curIdx, :) = curSynIds([curSynIdx, curNeighIdx]);
end

%% Generate skeleton for annotations
skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

for curIdx = 1:synPairCount
    curAxonId = randAxonIds(curIdx);
    curParentId = conn.axonMeta.parentId(curAxonId);
    curTreeIds = skel.numTrees() + (1:3);
    
    skel = Superagglos.toSkel(axons(curParentId), skel);
    skel.names{end} = sprintf('Axon %d', curAxonId);
    
    curSynIds = randSynIds(curIdx, :);
    curSynAgglos = cellfun( ...
        @vertcat, ...
        syn.synapses.presynId(curSynIds), ...
        syn.synapses.postsynId(curSynIds), ...
        'UniformOutput', false);
    curSynAgglos = cellfun( ...
        @(ids) segPoints(ids, :), ...
        curSynAgglos, 'UniformOutput', false);
    
    skel = Skeleton.fromMST( ...
        curSynAgglos, param.raw.voxelSize, skel);
    skel.names((end - 1):end) = arrayfun( ...
        @(id) sprintf('Synapse %d', id), ...
        curSynIds, 'UniformOutput', false);
    
   [skel, curGroupId] = skel.addGroup(sprintf( ...
        'Axon %d, Synapses %d and %d', ...
        curAxonId, curSynIds(1), curSynIds(2)));
    skel = skel.addTreesToGroup(curTreeIds, curGroupId);
end
