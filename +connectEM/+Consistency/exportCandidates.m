% This script exports axon-dendrite pairs that are coupled by exactly N
% synapses according to the connectome. These candidates are then inspected
% in webKNOSSOS and manually labelled as true or false.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outputDir = '/home/amotta/Desktop/most-distant';
dendriteFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_03_v2.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
asiRunId = '20190227T082543';

synCount = 2;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

points = Seg.Global.getSegToPointMap(param);

[conn, syn] = connectEM.Consistency.loadConnectome(param);

% Load axon-spine interfaces
[curDir, curAsiFile] = fileparts(connFile);
curAsiFile = sprintf('%s__%s_asiT.mat', curAsiFile, asiRunId);
curAsiFile = fullfile(curDir, curAsiFile);

asiT = load(curAsiFile);
asiT = asiT.asiT;

asiT = asiT(asiT.area > 0, :);
asiT = connectEM.Consistency.Calibration.apply(asiT);

% NOTE(amotta): Axons were subject to segment pick-up prior to connectome
% inference. The super-agglomerates are thus incomplete vis-à-vis the
% connectome.
%   The dendrite super-agglomerates were reduced to their segment-based
% components for the connectome. It's thus entirely valid to use their
% skeleton representation here. Let's do this (as it avoids the painfully
% slow minimum spanning tree calculation).
dendrites = load(dendriteFile);
dendrites = dendrites.dendrites(conn.denMeta.parentId);

%% Restrict to N-fold coupled neurites
clear cur*;

curAsiT = asiT;
curAsiT = curAsiT( ...
    curAsiT.type == 'PrimarySpine' & ismember( ...
    curAsiT.axonClass, {'Corticocortical', 'Thalamocortical'}), :);

pairT = curAsiT(:, {'preAggloId', 'postAggloId'});
[~, pairIds, pairSynCount] = unique(pairT, 'rows');

pairT = pairT(pairIds, :);
pairT.synCount = accumarray(pairSynCount, 1);
pairT(pairT.synCount ~= synCount, :) = [];
pairT.synCount = [];
clear pairSynCount;

[~, curAsiT.axonDendPairId] = ismember( ...
	curAsiT(:, {'preAggloId', 'postAggloId'}), pairT, 'rows');
curAsiT(~curAsiT.axonDendPairId, :) = [];

curAsiT = sortrows(curAsiT, 'axonDendPairId');
pairT.synIds = transpose(reshape(curAsiT.id, synCount, []));

%% Generate NML files
clear cur*;

rng(0);
randIds = randperm(height(pairT));
randIds = randIds(1:20);

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = Skeleton.setDescriptionFromRunInfo(skel, info);

for curIdx = 1:numel(randIds)
    curId = randIds(curIdx);
    curAxonId = pairT.preAggloId(curId);
    curDendId = pairT.postAggloId(curId);
    curSynIds = pairT.synIds(curId, :);
    
    curAxon = conn.axons{curAxonId};
    curDend = dendrites(curDendId);
    
    curSyns = syn.synapses(curSynIds, :);
    curSyns = cellfun( ...
        @(pre, post) unique(cat(1, pre, post)), ...
        curSyns.presynId, curSyns.postsynId, ...
        'UniformOutput', false);
    
    curSkel = skel;
    curSkel = Skeleton.fromMST( ...
        points(curAxon, :), param.raw.voxelSize, curSkel);
    curSkel.names{end} = sprintf('Axon %d', curAxonId);
    curSkel.colors{end} = [1, 0, 0, 1];
    
    curSkel = curSkel.addTree( ...
        sprintf('Dendrite %d', curDendId), ...
        curDend.nodes(:, 1:3), curDend.edges, ...
        [0, 0, 1, 1]);
    
    curSkel = Skeleton.fromMST( ...
        cellfun( ...
            @(ids) points(ids, :), curSyns, ...
            'UniformOutput', false), ...
        param.raw.voxelSize, curSkel);
    
    curMask = cat(1, false(2, 1), true(synCount, 1));
    curSkel.names(curMask) = arrayfun( ...
        @(id) sprintf('Synapse %d', id), ...
        curSynIds, 'UniformOutput', false);
    curSkel.colors(curMask) = {[0, 1, 0, 1]};
    
    curSkelFile = sprintf( ...
        '%0*d_axon-%d_dendrite_%d.nml', ...
        ceil(log10(1 + numel(randIds))), ...
        curIdx, curAxonId, curDendId);
    curSkel.write(fullfile(outputDir, curSkelFile));
end
