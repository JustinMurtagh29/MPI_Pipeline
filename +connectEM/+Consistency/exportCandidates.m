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

%% configuration
rootDir = '/home/amotta/Desktop'; % TODO(amotta)
outputDir = '/home/amotta/Desktop';

synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');

info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

param.saveFolder = rootDir; % TODO(amotta)
points = Seg.Global.getSegToPointMap(param);

syn = load(synFile);
conn = load(connFile);

%% limit synapses
% from `+connectEM/+Connectome/plotSynapseSizeCorrelations.m`
synT = table;
synT.id = cell2mat(conn.connectome.synIdx);
synT.area = cell2mat(conn.connectomeMeta.contactArea);
synT.isSpine = syn.isSpineSyn(synT.id);

synT.preAggloId = repelem( ...
    conn.connectome.edges(:, 1), ...
    cellfun(@numel, conn.connectome.synIdx));

synT.postAggloId = repelem( ...
    conn.connectome.edges(:, 2), ...
    cellfun(@numel, conn.connectome.synIdx));

% limit to spine synapses
synT(~synT.isSpine, :) = [];
synT.isSpine = [];

% remove duplicate entries
[~, uniRows, uniCount] = unique(synT.id);
synT = synT(uniRows, :);

% remove synapses occuring multiple times
% (i.e., between two or more pairs of neurites)
synT.occurences = accumarray(uniCount, 1);
synT(synT.occurences > 1, :) = [];
synT.occurences = [];

%% export multiply coupled neurite
[prePostPair, ~, prePostMult] = unique( ...
    synT(:, {'preAggloId', 'postAggloId'}), 'rows');
prePostMult = accumarray(prePostMult, 1);
prePostPair(prePostMult ~= 4, :) = [];

%% generate NML files
skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);

skelDesc = sprintf('%s (%s)', mfilename, info.git_repos{1}.hash);
skel = skel.setDescription(skelDesc);

rng(0);
randIds = randperm(size(prePostPair, 1), 25);

for curIdx = 1:numel(randIds)
    curId = randIds(curIdx);
    curAxonId = prePostPair.preAggloId(curId);
    curDendId = prePostPair.postAggloId(curId);
    
    curAxon = conn.axons{curAxonId};
    curDend = conn.dendrites{curDendId};
    
    curSkel = skel;
    curSkel = Skeleton.fromMST( ...
        points(curAxon, :), param.raw.voxelSize, curSkel);
    curSkel.names{end} = sprintf('Axon %d', curAxonId);
    curSkel.colors{end} = [1, 0, 0, 1];
    
    curSkel = Skeleton.fromMST( ...
        points(curDend, :), param.raw.voxelSize, curSkel);
    curSkel.names{end} = sprintf('Dendrite %d', curDendId);
    curSkel.colors{end} = [0, 0, 1, 1];
    
    curSkelFile = sprintf( ...
        '%0*d_axon-%d_dendrite_%d.nml', ...
        ceil(log10(1 + numel(randIds))), ...
        curIdx, curAxonId, curDendId);
    curSkel.write(fullfile(outputDir, curSkelFile));
end
