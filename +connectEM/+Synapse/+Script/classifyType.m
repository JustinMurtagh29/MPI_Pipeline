% Classifies synapses into
% * soma synapses
% * shaft synapses
% * primary spine synapse
% * secondary spine synapse
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_spine_attachment.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');

% As in Benedikt's +L4/updateParamsToNewestFiles.m
% Commit hash `590d8538d65463151d43491e2446e25ca11dd5f6`
graphFile = fullfile(rootDir, 'graphNew.mat');
synScoreFile = fullfile(rootDir, 'globalSynScores.mat');

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

% Prepare for load synapse scores
param.svg.graphFile = graphFile;
param.svg.synScoreFile = synScoreFile;

% Synapse scores
maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);
graph = Seg.IO.loadGraph(param, false);

% Spine heads
shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

conn = load(connFile);
syn = load(synFile);

%% Classify synapses
somaAgglos = conn.denMeta.targetClass == 'Somata';
somaAgglos = conn.dendrites(somaAgglos);

syn.synapses.type = ...
    connectEM.Synapse.classifyType( ...
        param, syn.synapses, graph.synScores, ...
        shAgglos, somaAgglos, conn.axons);

%% Build spine head table
shLUT = Agglo.buildLUT(maxSegId, shAgglos);

shSynT = syn.synapses;
shSynT.id = reshape( ...
    1:size(shSynT), [], 1);

[~, shSynT.typeId] = ismember( ...
    shSynT.type, {'PrimarySpine', 'SecondarySpine'});
shSynT(~shSynT.typeId, :) = [];

shSynT.shId = cellfun( ...
    @(ids) max(shLUT(ids)), ...
    shSynT.postsynId);
assert(all(shSynT.shId));

shSynIds = accumarray( ...
   [shSynT.shId, shSynT.typeId], shSynT.id, ...
   [numel(shAgglos), 2], @(ids) {ids}, {zeros(0, 1)});

%% Look at examples in webKNOSSOS
rng(0);

randIds = find(all(cellfun(@numel, shSynIds), 2));
randIds = randIds(randperm(numel(randIds)));
randIds = randIds(1:25);

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

for curIdx = 1:numel(randIds)
    curShId = randIds(curIdx);
    curSynIds = shSynIds(curShId, :);
    
    curSynTypes = repelem( ...
        [true; false], cellfun(@numel, curSynIds(:)));
    curSynIds = cell2mat(curSynIds(:));
    
    curSynAgglos = cellfun( ...
        @vertcat, ...
        syn.synapses.presynId(curSynIds), ...
        syn.synapses.postsynId(curSynIds), ...
        'UniformOutput', false);
    curShAgglo = shAgglos(curShId);
    
    curAgglos = cellfun( ...
        @(segIds) segPoints(segIds, :), ...
        [curShAgglo; curSynAgglos], ...
        'UniformOutput', false);
    
    curNames = arrayfun( ...
        @(id) sprintf('Synapse %d', id), ...
        curSynIds, 'UniformOutput', false);
    curNames(curSynTypes) = strcat( ...
        curNames(curSynTypes), {' (Primary)'});
    curNames(~curSynTypes) = strcat( ...
        curNames(~curSynTypes), {' (Secondary)'});
    
    curNames = [ ...
       {sprintf('Spine head %d', curShId)}; curNames(:)];
    
    curPrefix = sprintf('%0*d.', ceil(log10(1 + numel(randIds))), curIdx);
    curNames = strcat(curPrefix, {' '}, curNames);
    
    skel = Skeleton.fromMST(curAgglos, param.raw.voxelSize, skel);
    skel.names((end - numel(curNames) + 1):end) = curNames;
end
