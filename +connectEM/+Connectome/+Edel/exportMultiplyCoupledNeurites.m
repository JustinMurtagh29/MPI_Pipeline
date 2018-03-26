% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outputDir = '/home/amotta/Desktop';

synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v5_ax18a_deWC01wSp_withShId.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_ax18a_deWC01wSp_v4.mat');

edgeFile = fullfile(rootDir, 'globalEdges.mat');
edgeMapFile = fullfile(rootDir, 'SVGDB', 'agglos', 'ax18a_deWC01wSp', 'edges.mat');
segMapFile = fullfile(rootDir, 'SVGDB', 'agglos', 'ax18a_deWC01wSp', 'eClass.mat');

info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

points = Seg.Global.getSegToPointMap(param);

syn = load(synFile);
conn = load(connFile);

edges = load(edgeFile, 'edges');
edges = edges.edges;

edgeMap = load(edgeMapFile, 'borderMapping');
edgeMap = edgeMap.borderMapping;

segMap = load(segMapFile);
segMap = segMap.segMapping;

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

% limit to shaft synapses
synT(synT.isSpine, :) = [];
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
prePostPair(prePostMult ~= 2, :) = [];

%% generate NML files
skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);

skelDesc = sprintf('%s (%s)', mfilename, info.git_repos{1}.hash);
skel = skel.setDescription(skelDesc);

rng(0);
randIds = randperm(size(prePostPair, 1), 25);

for curIdx = 1% :numel(randIds)
    curId = randIds(curIdx);
    curAxonId = prePostPair.preAggloId(curId);
    curDendId = prePostPair.postAggloId(curId);
    
    curSynIds = synT.id( ...
        synT.preAggloId == curAxonId ...
      & synT.postAggloId == curDendId);
    curSynEdges = syn.synapses.edgeIdx(curSynIds);
    
    % Translate edges to SegEM-segmentation
    curSynSegIds = cellfun(@(ids) ...
        unique(edges(ismember(edgeMap, ids), :)), ...
        curSynEdges, 'UniformOutput', false);
    
    curAxon = find(ismember(segMap, conn.axons{curAxonId}));
    curDend = find(ismember(segMap, conn.dendrites{curDendId}));
    
    curSkel = skel;
    curSkel = Skeleton.fromMST( ...
        points(curAxon, :), param.raw.voxelSize, curSkel);
    curSkel.names{end} = sprintf('Axon %d', curAxonId);
    curSkel.colors{end} = [1, 0, 0, 1];
    
    curSkel = Skeleton.fromMST( ...
        points(curDend, :), param.raw.voxelSize, curSkel);
    curSkel.names{end} = sprintf('Dendrite %d', curDendId);
    curSkel.colors{end} = [0, 0, 1, 1];
    
    for curSynIdx = 1:numel(curSynIds)
        curSkel = Skeleton.fromMST( ...
            points(curSynSegIds{curSynIdx}, :), ...
            param.raw.voxelSize, curSkel);
        curSkel.names{end} = sprintf( ...
            'Synapse %d', curSynIds(curSynIdx));
        curSkel.colors{end} = [0, 1, 0, 1];
    end
    
    curSkelFile = sprintf( ...
        '%0*d_axon-%d_dendrite_%d.nml', ...
        ceil(log10(1 + numel(randIds))), ...
        curIdx, curAxonId, curDendId);
    curSkel.write(fullfile(outputDir, curSkelFile));
end
