% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');

outputDir = '';

exportConfigs = struct;
exportConfigs(1).name = 'large';
exportConfigs(1).cvRange = [0, 0.4];
exportConfigs(1).avgLogAsiRange = [-0.14, 0.23];

exportConfigs(2).name = 'small';
exportConfigs(2).cvRange = [0, 0.5];
exportConfigs(2).avgLogAsiRange = [-0.97, -0.6];

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

voxelSize = param.raw.voxelSize;
points = Seg.Global.getSegToPointMap(param);

[conn, syn, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

[~, synToSynFile] = fileparts(connFile);
synToSynFile = sprintf('%s_synToSynDists.mat', synToSynFile);
synToSynFile = fullfile(fileparts(connFile), synToSynFile);
synToSyn = load(synToSynFile);

% Loading spine head agglomerates
shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

% Loading augmented graph
graph = Graph.load(rootDir);
graph(~graph.borderIdx, :) = [];

borderAreas = fullfile(rootDir, 'globalBorder.mat');
borderAreas = load(borderAreas, 'borderArea2');
borderAreas = borderAreas.borderArea2;

graph.borderArea = borderAreas(graph.borderIdx);
clear borderAreas;

% HACK(amotta): For some reason there exist borders, for which
% `physicalBorderArea2` is zero. This seems wrong.
%   In order not to be affected by this issue, let's set the area of these
% borders to NaN. This will result in a total axon-spine interface area of
% NaN, which we can remove by brute force later on.
%
% Corresponding issue on GitLab:
% https://gitlab.mpcdf.mpg.de/connectomics/auxiliaryMethods/issues/16
graph.borderArea(~graph.borderArea) = nan;
graph(:, {'prob', 'borderIdx'}) = [];

%% Build axon-spine interface areas
asiT = ...
    connectEM.Connectome.buildAxonSpineInterfaces( ...
        param, graph, shAgglos, conn, syn);
asiT(asiT.type ~= 'PrimarySpine', :) = [];

% HACK(amotta): This is the counter-part to the above hack. (This is not
% really needed, since `buildAxonSpineInterface` does exactly the same
% thing internally. But let's be explicit and future-proof.)
asiT(isnan(asiT.area), :) = [];

%% Find synapse pairs
clear cur*;
curSynIdPairs = ...
    connectEM.Consistency.buildPairConfigs( ...
        asiT, struct('synIds', reshape(1:height(asiT), [], 1)));
curSynIdPairs = curSynIdPairs(1).synIdPairs;

pairT = asiT(curSynIdPairs(:, 1), {'preAggloId', 'postAggloId'});
pairT.synIds = asiT.id(curSynIdPairs);
pairT.areas = asiT.area(curSynIdPairs);
clear curSynIdPairs;

pairT.cv = std(pairT.areas, 0, 2) ./ mean(pairT.areas, 2);
pairT.avgLogAreas = mean(log10(pairT.areas), 2);

curInRange = @(r, v) r(1) <= v & v <= r(end);

for curIdx = 1:numel(exportConfigs)
    curConfig = exportConfigs(curIdx);
    
    curMask = ...
        curInRange(curConfig.cvRange, pairT.cv) ...
      & curInRange(curConfig.avgLogAsiRange, pairT.avgLogAreas);
    exportConfigs(curIdx).pairIds = reshape(find(curMask), [], 1);
end

%% Generate NML files
clear cur*;

if ~isempty(outputDir)
    skel = skeleton();
    skel = Skeleton.setParams4Pipeline(skel, param);

    skelDesc = sprintf('%s (%s)', mfilename, info.git_repos{1}.hash);
    skel = skel.setDescription(skelDesc);

    for curConfig = exportConfigs
        rng(0);
        randIds = curConfig.pairIds;
        randIds = randIds(randperm(numel(randIds)));
        randIds = randIds(1:50);

        curOutputDir = fullfile(outputDir, curConfig.name);
        if ~exist(curOutputDir, 'dir'); mkdir(curOutputDir); end

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
                points(curAxon, :), voxelSize, curSkel);
            curSkel.names{end} = sprintf('Axon %d', curAxonId);
            curSkel.colors{end} = [1, 0, 0, 1];

            curSkel = Skeleton.fromMST( ...
                    points(curDend, :), voxelSize, curSkel);
            curSkel.names{end} = sprintf('Dendrite %d', curDendId);
            curSkel.colors{end} = [0, 0, 1, 1];

            curSkel = Skeleton.fromMST( ...
                cellfun( ...
                    @(ids) points(ids, :), curSyns, ...
                    'UniformOutput', false), ...
                param.raw.voxelSize, curSkel);

            curMask = cat(1, false(2, 1), true(numel(curSynIds), 1));
            curSkel.names(curMask) = arrayfun( ...
                @(id) sprintf('Synapse %d', id), ...
                curSynIds, 'UniformOutput', false);
            curSkel.colors(curMask) = {[0, 1, 0, 1]};

            curSkelFile = sprintf( ...
                '%0*d_axon-%d_dendrite_%d.nml', ...
                ceil(log10(1 + numel(randIds))), ...
                curIdx, curAxonId, curDendId);
            curSkel.write(fullfile(curOutputDir, curSkelFile));
        end
    end
end

%% Look at potential distance-dependent effects
clear cur*;
for curConfig = exportConfigs
    curSynT = asiT;
    curSynToSyn = synToSyn;
    
    curPairConfig = struct;
    curPairConfig.title = curConfig.name;
   [~, curPairConfig.synIdPairs] = ismember( ...
       pairT.synIds(curConfig.pairIds, :), asiT.id);
    
    connectEM.Consistency.plotVariabilityVsDistance( ...
        curSynT, curSynToSyn, curPairConfig, 'maxDistUm', 50);
end
