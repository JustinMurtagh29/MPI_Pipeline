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
outputDir = '/home/amotta/Desktop';

synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_b_linearized_ax_spine_syn_clust.mat');

synCount = 2;
synType = 'spine';

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

points = Seg.Global.getSegToPointMap(param);

syn = load(synFile);
conn = load(connFile);

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

[axonDendPair, ~, pairSynCount] = unique( ...
    synT(:, {'preAggloId', 'postAggloId'}), 'rows');

pairSynCount = accumarray(pairSynCount, 1);
axonDendPair(pairSynCount ~= synCount, :) = [];
clear pairSynCount;

%% Generate NML files
skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);

skelDesc = sprintf('%s (%s)', mfilename, info.git_repos{1}.hash);
skel = skel.setDescription(skelDesc);

openIds = 1:size(axonDendPair, 1);

% First round
rng(0);
randIds = openIds(randperm(numel(openIds), 20));
openIds = setdiff(openIds, randIds);

% Second round
rng(0);
randIds = openIds(randperm(numel(openIds), 20));

for curIdx = 1:numel(randIds)
    curId = randIds(curIdx);
    curAxonId = axonDendPair.preAggloId(curId);
    curDendId = axonDendPair.postAggloId(curId);
    
    curSynIds = synT.id( ...
        synT.preAggloId == curAxonId ...
      & synT.postAggloId == curDendId);
    
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
