% This script exports axon-dendrite pairs that are coupled by exactly N
% synapses according to the connectome. These candidates are then inspected
% in webKNOSSOS and manually labelled as true or false.
%
% In the "true" case, the pre- and postsynaptic volumes are manually
% reconstructed by picking up the corresponding segments. This will give a
% good estimate of the true contact area.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outputDir = '/tmpscratch/amotta/l4/2018-06-07-distance-dependence-of-synaptic-consistency-check/closer-than-5um';

connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

synCount = 2;
synType = 'spine';

distType = 'preDist';
minDist = [];
maxDist = 5E3;

info = Util.runInfo();

%% Check and complete configuration
calcDist = false;
if ~isempty(minDist) || ~isempty(maxDist)
    % Distance calculation only works for synapse pairs
    assert(ismember(distType, {'preDist', 'postDist'}));
    assert(synCount == 2);
    calcDist = true;
end

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

points = Seg.Global.getSegToPointMap(param);
[conn, syn] = connectEM.Connectome.load(param, connFile);

if calcDist
   [~, synToSynFile] = fileparts(connFile);
    synToSynFile = sprintf('%s_synToSynDists.mat', synToSynFile);
    synToSynFile = fullfile(fileparts(connFile), synToSynFile);
    synToSyn = load(synToSynFile);
end

%% Restrict to N-fold coupled neurites
synT = connectEM.Connectome.buildSynapseTable(conn, syn);

switch lower(synType)
    case 'spine'
        synT(~synT.isSpine, :) = [];
    case 'shaft'
        targetClass = conn.denMeta.targetClass(synT.postAggloId);
        synT(synT.isSpine | (targetClass == 'Somata'), :) = [];
        clear targetClass;
    otherwise
        error('Unknown type "%s"', synType)
end

pairT = synT(:, {'preAggloId', 'postAggloId'});
[~, pairIds, pairSynCount] = unique(pairT, 'rows');

pairT = pairT(pairIds, :);
pairT.synCount = accumarray(pairSynCount, 1);
pairT(pairT.synCount ~= synCount, :) = [];
pairT.synCount = [];
clear pairSynCount;

[~, synT.axonDendPairId] = ismember( ...
	synT(:, {'preAggloId', 'postAggloId'}), pairT, 'rows');
synT(~synT.axonDendPairId, :) = [];

synT = sortrows(synT, 'axonDendPairId');
pairT.synIds = transpose(reshape(synT.id, synCount, []));

%% Calculate pair-wise distances
if calcDist
    pairT.preDist = zeros(size(pairT.preAggloId));
    pairT.postDist = zeros(size(pairT.postAggloId));
    
    for curIdx = 1:height(pairT)
        curPreAggloId = pairT.preAggloId(curIdx);
        curPostAggloId = pairT.postAggloId(curIdx);
        curSynIds = pairT.synIds(curIdx, :);

        % Distance along axonal side
        curPreSynIds = synToSyn.axonSynIds{curPreAggloId};
       [~, curPreSynIds] = ismember(curSynIds, curPreSynIds);
        curPreDist = synToSyn.axonSynToSynDists{curPreAggloId};
        curPreDist = curPreDist(curPreSynIds(1), curPreSynIds(2));

        % Distance along axonal side
        curPostSynIds = synToSyn.dendSynIds{curPostAggloId};
       [~, curPostSynIds] = ismember(curSynIds, curPostSynIds);
        curPostDist = synToSyn.dendSynToSynDists{curPostAggloId};
        curPostDist = curPostDist(curPostSynIds(1), curPostSynIds(2));

        pairT.preDist(curIdx) = curPreDist;
        pairT.postDist(curIdx) = curPostDist;
    end
    
    if ~isempty(maxDist); pairT(pairT.(distType) > maxDist, :) = []; end
    if ~isempty(minDist); pairT(pairT.(distType) < minDist, :) = []; end
end

%% Generate NML files
skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);

skelDesc = sprintf('%s (%s)', mfilename, info.git_repos{1}.hash);
skel = skel.setDescription(skelDesc);

rng(0);
randIds = randperm(height(pairT));
randIds = randIds(1:25);

for curIdx = 1:numel(randIds)
    curId = randIds(curIdx);
    curAxonId = pairT.preAggloId(curId);
    curDendId = pairT.postAggloId(curId);
    curSynIds = pairT.synIds(curId, :);
    
    curAxon = conn.axons{curAxonId};
    curDend = conn.dendrites{curDendId};
    
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
    
    curSkel = Skeleton.fromMST( ...
        points(curDend, :), param.raw.voxelSize, curSkel);
    curSkel.names{end} = sprintf('Dendrite %d', curDendId);
    curSkel.colors{end} = [0, 0, 1, 1];
    
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
