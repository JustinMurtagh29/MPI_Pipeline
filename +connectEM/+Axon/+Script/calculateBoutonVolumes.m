% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);
segPoints = param.raw.voxelSize .* segPoints;

segVols = Seg.Global.getSegToSizeMap(param);
segVols = prod(param.raw.voxelSize) .* segVols;

[conn, syn] = connectEM.Connectome.load(param, connFile);

[connDir, connName] = fileparts(connFile);
interSynFile = sprintf('%s_intersynapse_v2.mat', connName);
interSynFile = fullfile(connDir, interSynFile);
interSyn = load(interSynFile);

%% Testing
synIds = connectEM.Axon.getSynapses(conn, syn);
boutonIds = connectEM.Axon.clusterSynapsesIntoBoutons(synIds, interSyn);

%%
connectEM.Axon.buildBoutonAgglos( ...
    segPoints, conn, syn, interSyn, ...
    'showProgressBar', true);

%% For axon
clear cur*;
curDistThresh = 1000;

curAxonId = 4;
curSegIds = conn.axons{curAxonId};
curKdTree = KDTreeSearcher(segPoints(curSegIds, :));

curGraph = SuperAgglo.fromAgglo({curSegIds}, segPoints, 'mst');
curGraph = graph( ...
    curGraph.edges(:, 1), curGraph.edges(:, 2), ...
    sqrt(sum(( ...
        curGraph.nodes(curGraph.edges(:, 1), 1:3) ...
      - curGraph.nodes(curGraph.edges(:, 2), 1:3)) .^ 2, 2)), ...
	size(curGraph.nodes, 1));

curSynIds = synIds{curAxonId};
curBoutonIds = boutonIds{curAxonId};
curBoutonCount = max(curBoutonIds);

for curBoutonId = 1:curBoutonCount
    curSynSegIds = curSynIds(curBoutonIds == curBoutonId);
    curSynSegIds = syn.synapses.presynId(curSynSegIds);
    curSynSegIds = unique(cell2mat(curSynSegIds));
    
    curBoutonSegIds = rangesearch( ...
        curKdTree, segPoints(curSynSegIds, :), curDistThresh);
    curBoutonSegIds = cell2mat(reshape(curBoutonSegIds, 1, []));
    curBoutonSegIds = unique(curBoutonSegIds(:));
    
   [~, curSynSegIds] = ismember(curSynSegIds, curSegIds);
    curSynSegIds = setdiff(curSynSegIds, 0);
    
    curDists = distances(curGraph, curSynSegIds, curBoutonSegIds);
    curBoutonSegIds = curBoutonSegIds(any(curDists < curDistThresh, 1));
    curBoutonVol = sum(segVols(curBoutonSegIds));
end
