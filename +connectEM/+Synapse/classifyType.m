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

% Somata
conn = load(connFile);

somaAgglos = (conn.denMeta.targetClass == 'Somata');
somaAgglos = conn.dendrites(somaAgglos);
axonAgglos = conn.axons;

% Synapses
syn = load(synFile);
synapses = syn.synapses;

%%
% Building look-up tables
% Soma class dominates over spine heads
axonLUT = Agglo.buildLUT(maxSegId, axonAgglos);
somaLUT = Agglo.buildLUT(maxSegId, somaAgglos);
shLUT = Agglo.buildLUT(maxSegId, shAgglos);
shLUT(somaLUT ~= 0) = 0;

synapses.id = reshape( ...
    1:size(synapses, 1), [], 1);
synapses.synScores = cellfun( ...
    @(ids) max(graph.synScores(ids, :), [], 2), ...
    synapses.edgeIdx, 'UniformOutput', false);

% Sanity check
% Synapse threshold was -1.67
assert(all(cellfun(@min, synapses.synScores) > -1.67));
synapses.maxSynScore = cellfun(@max, synapses.synScores);

synapses.shId = cellfun( ...
    @(segIds) max(shLUT(segIds)), ...
    synapses.postsynId);
synapses.axonId = cellfun( ...
    @(segIds) max(axonLUT(segIds)), ...
    synapses.presynId);

%% Separate primary from secondary spine synapses
shSynIds = accumarray( ...
    1 + synapses.shId, synapses.id, ...
   [1 + numel(shAgglos), 1], @(ids) {ids}, {zeros(0, 1)});
shSynIds = shSynIds(2:end);

% Sort synapse by SynEM scores
[~, sortedIds] = cellfun( ...
    @(synIds) sort(synapses.maxSynScore(synIds)), ...
    shSynIds, 'UniformOutput', false);
shSynIds = cellfun( ...
    @(synIds, sortIds) synIds(sortIds), ...
    shSynIds, sortedIds, 'UniformOutput', false);
clear sortedIds;

shT = table;
shT.id = reshape( ...
    1:numel(shAgglos), [], 1);

shT(cellfun(@isempty, shSynIds), :) = [];
shSynIds(cellfun(@isempty, shSynIds)) = [];

shT.priSynId = cellfun( ...
    @(synIds) synIds(end), shSynIds);
shT.secSynIds = cellfun( ...
    @(synIds) reshape(synIds(1:(end - 1)), [], 1), ...
    shSynIds, 'UniformOutput', false);

%% Handle split spine synapses
priSpineSynIds = ismember(synapses.id, shT.priSynId);
priSpineSynIds = synapses(priSpineSynIds, :);

priSpineSynIds = ismember( ...
    synapses(:,  {'shId', 'axonId'}), ...
    priSpineSynIds(:, {'shId', 'axonId'}), 'rows');
priSpineSynIds = synapses.id(priSpineSynIds);

secSpineSynIds = cell2mat(shT.secSynIds);
secSpineSynIds = setdiff(secSpineSynIds, priSpineSynIds);

somaSynIds = cellfun( ...
    @(segIds) max(somaLUT(segIds)), ...
    synapses.postsynId);
somaSynIds = find(somaSynIds);

% Sanity check
assert(isempty(intersect(priSpineSynIds, secSpineSynIds)));

%%
synapses.type(:) = {'Shaft'};
synapses.type(priSpineSynIds) = {'PrimarySpine'};
synapses.type(secSpineSynIds) = {'SecondarySpine'};
synapses.type(somaSynIds) = {'Soma'};
synapses.type = categorical(synapses.type);

%%


%% Look at examples in webKNOSSOS
%{
rng(0);

randIds = find(~cellfun(@isempty, shT.secSynIds));
randIds = randIds(randperm(numel(randIds)));
randIds = randIds(1:25);

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

for curIdx = 1:numel(randIds)
    curShT = shT(randIds(curIdx), :);
    curShAgglo = shAgglos(curShT.id);
    
    curSynIds = [ ...
        curShT.priSynId; ...
        curShT.secSynIds{1}];
    curSynAgglos = cellfun( ...
        @vertcat, ...
        synapses.presynId(curSynIds), ...
        synapses.postsynId(curSynIds), ...
        'UniformOutput', false);
    
    curAgglos = cellfun( ...
        @(segIds) segPoints(segIds, :), ...
        [curShAgglo; curSynAgglos], ...
        'UniformOutput', false);
    
    curNames = [ ...
       {sprintf( ...
            'Spine head %d', curShT.id)}; ...
        arrayfun( ...
            @(id) sprintf('Synapse %d', id), ...
            curSynIds, 'UniformOutput', false)];
	curNames(2) = strcat(curNames(2), {' (Primary)'});
    curNames(3:end) = strcat(curNames(3:end), {' (Secondary)'});
    
    curPrefix = sprintf('%0*d.', ceil(log10(1 + numel(randIds))), curIdx);
    curNames = strcat(curPrefix, {' '}, curNames);
    
    skel = Skeleton.fromMST(curAgglos, param.raw.voxelSize, skel);
    skel.names((end - numel(curNames) + 1):end) = curNames;
end
%}
