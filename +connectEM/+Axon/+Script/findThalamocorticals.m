% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a.mat');
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3.mat');

%% loading data
conn = load(connFile);
syn = load(synFile);

% load dataset parameters
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

% load segment positions
points = Seg.Global.getSegToPointMap(param);

% load border indices for edges
% TODO(amotta): make sure that Benedikt used the same graph!
segGraph = fullfile(rootDir, 'graph.mat');
segGraph = load(segGraph, 'borderIdx');

% load size and c. o. m. of borders
borders = fullfile(rootDir, 'globalBorder.mat');
borders = load(borders, 'borderCoM', 'borderArea2');

%% calculate synapse positions
synapses = syn.synapses;
synapses.ontoSpine = syn.isSpineSyn;

% convert edge indices to border indices
synapses.borderIdx = cellfun( ...
    @(ids) setdiff(segGraph.borderIdx(ids), 0), ...
    synapses.edgeIdx, 'UniformOutput', false);

synapses.pos = cell2mat(cellfun( ...
    @(i) ceil( ...
        sum( ...
            borders.borderArea2(i) ...
            .* double(borders.borderCoM(i, :)), 1 ...
        ) ./ sum(borders.borderArea2(i))), ...
    synapses.borderIdx, 'UniformOutput', false));

%% testing
% show an axon agglomerate an all its synapses
axonId = find(conn.axonMeta.synCount >= 10);

rng(0);
axonId = axonId(randperm(numel(axonId), 1));
axonSegIds = conn.axons{axonId};

% find synapses for axon
synapseIds = conn.connectome.edges(:, 1) == axonId;
synapseIds = conn.connectome.synIdx(synapseIds);
synapseIds = cell2mat(synapseIds);

% map synapses to segments
synapseNodeDist = pdist2( ...
    points(axonSegIds, :) .* param.raw.voxelSize, ...
    synapses.pos(synapseIds, :) .* param.raw.voxelSize);
[~, synapseNodeIds] = min(synapseNodeDist, [], 1);

% build MST representation
mstGraph = graphminspantree(sparse(squareform( ...
    pdist(points(axonSegIds, :) .* param.raw.voxelSize))));
mstGraph = graph(mstGraph, 'lower');

% calculate synapse-to-synapse distance (along MST)
[uniNodeIds, ~, synToUniNodes] = unique(synapseNodeIds);
synToSynDist = distances(mstGraph, uniNodeIds, uniNodeIds);
synToSynDist = synToSynDist(synToUniNodes, synToUniNodes);

% set self-distance to infinity
synToSynDist(1:(size(synToSynDist, 1) + 1):end) = inf;
synToSynDist(~synToSynDist) = eps;
closestSynDist = min(synToSynDist, [], 1);

figure;
histogram(closestSynDist / 1E3, linspace(0, 30, 31));

%{


skel = Skeleton.fromMST(points(axonSegIds, :), param.raw.voxelSize);
skel = skel.addNodesAsTrees(synapses.pos(synapseIds, :));
skel = Skeleton.setParams4Pipeline(skel, param);
skel.write('/home/amotta/Desktop/debug-synapse-locations.nml');
%}