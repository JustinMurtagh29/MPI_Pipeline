% Calculates the pairwise distance between synapses along the axon.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

[interSynFile, interSynName] = fileparts(connFile);
interSynName = sprintf('%s_intersynapse.mat', interSynName);
interSynFile = fullfile(interSynFile, interSynName);
clear interSynName;

info = Util.runInfo();

%% loading data
conn = load(connFile);
syn = load(conn.info.param.synFile);

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
out = struct;
out.info = info;
out.axonIds = find(conn.axonMeta.synCount >= 2);
out.axonPathLens = nan(size(out.axonIds));
out.synToSynDists = cell(size(out.axonIds));
out.synIds = cell(size(out.axonIds));

tic;
for idx = 1:numel(out.axonIds)
    axonId = out.axonIds(idx);
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
    out.axonPathLens(idx) = axonPathLen;
    out.synToSynDists{idx} = synToSynDist;
    out.synIds{idx} = synapseIds(:);
    
    Util.progressBar(idx, numel(out.axonIds));
end

% save results
Util.saveStruct(interSynFile, out);
