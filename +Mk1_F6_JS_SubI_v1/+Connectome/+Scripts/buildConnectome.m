%% script to generate the connectome for a segmentation v7
% See also: https://mhlablog.net/2018/07/05/l4-synapse-agglomeration
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
%         Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

axonFile = fullfile(param.saveFolder, 'aggloMat', '20191227T134319_ga_20191224T235355optimParams_agglomeration/20191227T220548_results.mat');
dendriteFile = fullfile(param.saveFolder, 'aggloMat', '20191227T134319_ga_20191224T235355optimParams_agglomeration/20191227T220548_results_auto-spines_v3.mat');

[~, curAxonName] = fileparts(axonFile);
[~, curDendName] = fileparts(dendriteFile);

synFile = fullfile( ...
    rootDir, 'connectome', sprintf( ...
        'SynapseAgglomerates__%s__%s__v1.mat', ...
        curAxonName, curDendName));
clear cur*;

info = Util.runInfo();
Util.showRunInfo(info);

%% Complete configuration
p = load(fullfile(rootDir, 'allParameter.mat'));
p = p.p;

p.svg = struct;
p.svg.edgeFile = fullfile(rootDir, 'globalEdges.mat');
p.svg.synScoreFile = fullfile(rootDir, 'globalSynapseScores.mat');
p.svg.segmentMetaFile = fullfile(rootDir, 'segmentMeta.mat');
p.svg.borderMetaFile = fullfile(rootDir, 'globalBorder.mat');
p.svg.graphFile = fullfile(rootDir, 'graph.mat');
p.svg.correspondenceFile = fullfile(rootDir, 'correspondences.mat');
p.svg.aggloPredFile = fullfile(rootDir, 'segmentAggloPredictions.mat');
p.svg.heuristicFile = fullfile(rootDir, 'heuristicResultComplete.mat');

p.connectome = struct;
p.connectome.saveFolder = fullfile(rootDir, 'connectome');

p.agglo = struct;
p.agglo.axonAggloFile = axonFile;
p.agglo.spineHeadFile = dendriteFile;

[~, axName] = fileparts(axonFile);
axName = strrep(axName, '_', '-');

[~, dendName] = fileparts(dendriteFile);
dendName = strrep(dendName, '_', '-');

[~, synName] = fileparts(synFile);
synName = strrep(synName, '_', '-');

outFile = fullfile(p.connectome.saveFolder, sprintf( ...
    'Connectome_%s_%s_%s.mat', axName, dendName, synName));
outFileSmall = fullfile(p.connectome.saveFolder, sprintf( ...
    'Connectome_%s_%s_%s_small.mat', axName, dendName, synName));

%% Loading agglos & synapses
% axons
Util.log('Loading axons.');
p.agglo.axonAggloFile = axonFile;
axons = Util.load(p.agglo.axonAggloFile, 'axons');
axons = arrayfun(@(a) double(a.segIds), axons, 'UniformOutput', false);
axonParentIds = reshape(1:numel(axons), size(axons));
assert(~any(cellfun(@isempty, axons)));

% dendrites
Util.log('Loading dendrites and target classes.');
dendrites = Util.load(dendriteFile, 'dendrites');
dendrites = arrayfun(@(d) d.nodes(:, 4), dendrites, 'UniformOutput', false);

% NOTE(amotta): We don't have any dendrite types yet...
% Let's consider them all to be "other dendrites" for now...
targetClass = num2cell(size(dendrites));
targetClass = repelem(categorical({'OtherDendrite'}), targetClass{:});

dendriteParentIds = find(targetClass ~= 'Ignore');
dendrites = dendrites(dendriteParentIds);
targetClass = targetClass(dendriteParentIds);

% target classes
idxSoma = targetClass == 'Somata';
idxWC = targetClass == 'WholeCell';
idxAD = targetClass == 'ApicalDendrite';
idxSD = targetClass == 'SmoothDendrite';
idxAIS = targetClass == 'AxonInitialSegment';
idxOther = targetClass == 'OtherDendrite';

% synapses
Util.log('Loading synapses.');
mSyn = load(synFile);
synapses = mSyn.synapses;

%% calculate full connectome

[ syn2Agglo, pre2syn, post2syn, connectome ] = ...
    L4.Synapses.synapsesAggloMapping( synapses, axons, dendrites );
connectomeFull = L4.Connectome.pairwise2Full(connectome, ...
    [length(axons), length(dendrites)]);


%% connectome synapse information

Util.log('Adding synapse meta information.');

% load border meta (must be wrt to edges in graph)
[graph, segmentMeta, borderMeta] = Seg.IO.loadGraph(p, false);
dend = load(p.svg.borderMetaFile, 'borderArea');
borderMeta.borderArea = nan(size(borderMeta.borderSize));
borderMeta.borderArea(~isnan(borderMeta.borderSize)) = dend.borderArea;

connectomeMeta = L4.Connectome.connectomeBorderMeta(connectome, ...
    synapses, borderMeta);


%% generate contactome

Util.log('Calculating contactome.')
contactome = L4.Agglo.findEdgesBetweenAgglos2(axons, dendrites, ...
    graph.edges);
contactomeMeta = L4.Connectome.contactomeBorderMeta(contactome, ...
    borderMeta);


%% target class connectome

Util.log('Target class mapping.');
denClasses = {'soma', 'whole cells', 'apical dendrites', ...
    'smooth dendrite', 'AIS', 'other'};
classConnectome = zeros(length(axons), length(denClasses));
classConnectome(:,1) = sum(connectomeFull(:, idxSoma), 2);
classConnectome(:,2) = sum(connectomeFull(:, idxWC), 2);
classConnectome(:,3) = sum(connectomeFull(:, idxAD), 2);
classConnectome(:,4) = sum(connectomeFull(:, idxSD), 2);
classConnectome(:,5) = sum(connectomeFull(:, idxAIS), 2);
classConnectome(:,6) = sum(connectomeFull(:, idxOther), 2);


%% axon statistics
% taken from connectEM.Connectome.augmentForMoritz

Util.log('Adding axon meta information.');
axonMeta = table;
axonMeta.id = reshape(1:numel(axons), [], 1);
axonMeta.parentId = axonParentIds;
axonMeta.synCount = accumarray( ...
    connectome.edges(:, 1), ...
    cellfun(@numel, connectome.synIdx), ...
   [numel(axons), 1], @sum, 0);
axonMeta.spineSynCount = accumarray( ...
    connectome.edges(:, 1), cellfun(@(ids) ...
    sum(mSyn.isSpineSyn(ids)), connectome.synIdx), ...
   [numel(axons), 1], @sum, 0);


%% dendrite statistics
% added by amotta

Util.log('Adding dendrite meta information.');
denMeta = table;
denMeta.id = reshape(1:numel(dendrites), [], 1);
denMeta.parentId = dendriteParentIds;
denMeta.targetClass = reshape(targetClass, [], 1);
denMeta.synCount = accumarray( ...
    connectome.edges(:, 2), ...
    cellfun(@numel, connectome.synIdx), ...
   [numel(dendrites), 1], @sum, 0);
denMeta.spineSynCount = accumarray( ...
    connectome.edges(:, 2), cellfun(@(ids) ...
    sum(mSyn.isSpineSyn(ids)), connectome.synIdx), ...
   [numel(dendrites), 1], @sum, 0);


%% save results

save(outFile, 'connectome', 'contactome', 'connectomeMeta', ...
    'contactomeMeta', 'axons', 'dendrites', 'classConnectome', ...
    'denClasses', 'axonMeta', 'denMeta', 'info');
Util.protect(outFile);

save(outFileSmall, 'connectome', 'connectomeMeta', 'classConnectome', ...
    'denClasses', 'axonMeta', 'denMeta', 'info');
Util.protect(outFileSmall);
