% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a.mat');
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3.mat');

[interSynFile, interSynName] = fileparts(connFile);
interSynName = sprintf('%s_intersynapse.mat', interSynName);
interSynFile = fullfile(interSynFile, interSynName);
clear interSynName;

info = Util.runInfo();

%% loading data
conn = load(connFile);
syn = load(synFile);

% load dataset parameters
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

% load segment positions
points = Seg.Global.getSegToPointMap(param);

%%
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

% calculate synapse locations
synapses.pos = cell2mat(cellfun( ...
    @(i) ceil( ...
        sum( ...
            borders.borderArea2(i) ...
            .* double(borders.borderCoM(i, :)), 1 ...
        ) ./ sum(borders.borderArea2(i))), ...
    synapses.borderIdx, 'UniformOutput', false));

%% calculate synapse-to-synapse distances
% only consider axons with at least two output synapses
axonIds = find(conn.axonMeta.synCount >= 2);
axonPathLens = nan(size(axonIds));
synToSynDists = cell(size(axonIds));

tic;
for idx = 1:numel(axonIds)
    axonId = axonIds(idx);
    axonSegIds = conn.axons{axonId};

    % find synapses for axon
    synapseIds = conn.connectome.edges(:, 1) == axonId;
    synapseIds = conn.connectome.synIdx(synapseIds);
    synapseIds = cell2mat(synapseIds);

    if numel(axonSegIds) > 1
        % map synapses to segments
        synapseNodeDist = pdist2( ...
            points(axonSegIds, :) .* param.raw.voxelSize, ...
            synapses.pos(synapseIds, :) .* param.raw.voxelSize);
       [~, synapseNodeIds] = min(synapseNodeDist, [], 1);
       
        % build MST representation of axon
        mstGraph = graphminspantree(sparse(squareform( ...
            pdist(points(axonSegIds, :) .* param.raw.voxelSize))));
        
        axonPathLen = sum(mstGraph(:));
        mstGraph = graph(mstGraph, 'lower');

        % calculate synapse-to-synapse distance (along MST)
       [uniNodeIds, ~, synToUniNodes] = unique(synapseNodeIds);
        synToSynDist = distances(mstGraph, uniNodeIds, uniNodeIds);
        synToSynDist = synToSynDist(synToUniNodes, synToUniNodes);
    else
        % It's possible for axon to consist of a single segment and still
        % to have multiple synapses. In this case the adjacency matrix
        % produced by `pdist` will be empty.
        axonPathLen = 0;
        synToSynDist = zeros(numel(synapseIds));
    end

    % set self-distance to infinity
    synToSynDist(1:(size(synToSynDist, 1) + 1):end) = inf;
    
    % store output
    axonPathLens(idx) = axonPathLen;
    synToSynDists{idx} = synToSynDist;
    
    Util.progressBar(idx, numel(axonIds));
end

% save results
Util.save(interSynFile, info, axonIds, axonPathLens, synToSynDists);

%% calculate inter-synapse distances
% Previously we've calculate the synapse-to-synapse distances along the
% axon. To calculate the inter-synapse distances we now need to construct a
% tree-like representation of the axon. Let's do this by calculating the
% minimal spanning-tree over the synapses.
data = load(interSynFile, 'axonIds', 'axonPathLens', 'synToSynDists');

conn.axonMeta.pathLen = nan(size(conn.axonMeta.id));
conn.axonMeta.interSynDists = cell(size(conn.axonMeta.id));

for curIdx = 1:numel(data.axonIds)
    curAxonId = data.axonIds(curIdx);
    curPathLen = data.axonPathLens(curIdx) / 1E3;
    
    curSynToSynDists = data.synToSynDists{curIdx};
    curSynToSynDists = curSynToSynDists ./ 1E3;
    
    % NOTE(amotta): Zeros in the adjacency matrix are interpreted as
    % missing edge. This is a problem since synapse-to-synapse distance
    % zero occurs naturally when multiple synapses originate from the same
    % segment. Let's instead set these zeros to the smallest possible
    % non-zero value.
    curSynToSynDists(~curSynToSynDists) = eps;
    
    % claculate inter-synapse distances
    curInterSynDists = graphminspantree(sparse(curSynToSynDists));
    curInterSynDists = nonzeros(curInterSynDists);
    
    conn.axonMeta.pathLen(curAxonId) = curPathLen;
    conn.axonMeta.interSynDists{curAxonId} = curInterSynDists;
end

%% orthogonal approach
minSynCount = 10;

% empirically found
minSynDensity = 0.25;
minSpineSynFrac = 0.5;
maxMedianInterSynDist = 2.2;

% spine synapse fraction
conn.axonMeta.spineFrac = ...
    conn.axonMeta.spineSynCount ...
    ./ conn.axonMeta.synCount;

% synapse density
conn.axonMeta.synDensity = ...
    conn.axonMeta.synCount ...
    ./ conn.axonMeta.pathLen;

% median inter-synapse distance
conn.axonMeta.medianInterSynDist = ...
    cellfun(@median, conn.axonMeta.interSynDists);

candMask = ...
    (conn.axonMeta.synCount >= minSynCount) ...
  & (conn.axonMeta.spineFrac >= minSpineSynFrac) ...
  & (conn.axonMeta.synDensity >= minSynDensity) ...
  & (conn.axonMeta.medianInterSynDist <= maxMedianInterSynDist);

% select random subset
rng(0);
candAxonIds = conn.axonMeta.id(candMask);
candAxonIds = candAxonIds(randperm(sum(candMask), 100));

skel = Skeleton.fromMST(cellfun( ...
    @(ids) points(ids, :), conn.axons(candAxonIds), ...
    'UniformOutput', false), param.raw.voxelSize);
skel.names = arrayfun( ...
    @(i, n) sprintf('%0*d. Axon %d', ...
        ceil(log10(1 + numel(candAxonIds))), i, n), ...
	reshape(1:numel(candAxonIds), [], 1), candAxonIds, ...
    'UniformOutput', false);

skel = Skeleton.setParams4Pipeline(skel, param);
skel.write('/home/amotta/Desktop/tc-axon-candidates.nml');