% Modified from +L4/+Synapses/+Scripts/synapseDetection_v3.m
% by Alessandro Motta <alessandro.motta@brain.mpg.de>
% 
% NOTE Indices in synapses are w.r.t. to the edges that were used during
%      synaptic interface clustering. This is done with the edges from the
%      graph, i.e. including correspondences.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
% Modified by: Sahil Loomba <sahil.loomba@brain.mpg.de>

% Use SVM predictions to overlay with SynEM predictions and get "better" synapse agglomerates

addpath(genpath('/u/sahilloo/repos/Benedikt/'));
rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
p = load(fullfile(rootDir, 'allParameter.mat'));

axonAggloFile = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/aggloMat/20191227T134319_ga_20191224T235355optimParams_agglomeration/20191227T220548_results.mat';
spineHeadFile = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/aggloMat/20191227T134319_ga_20191224T235355optimParams_agglomeration/20191227T220548_results_auto-spines_v3.mat';

[~, curAxonName] = fileparts(axonAggloFile);
[~, curSpineName] = fileparts(spineHeadFile);

outputFile = fullfile( ...
    rootDir, sprintf( ...
        'SynapseAgglomerates__%s__%s__v1.mat', ...
        curAxonName, curSpineName));
clear cur*;

info = Util.runInfo();
Util.showRunInfo(info);

%% load/set agglo files for heuristics (for reproducability)
p = p.p;

p.svg = struct;
p.svg.edgeFile = fullfile(rootDir, 'globalEdges.mat');
p.svg.synScoreFile = fullfile(rootDir, 'globalSynapseScores.mat');
p.svg.segmentMetaFile = fullfile(rootDir, 'segmentMeta.mat');
p.svg.borderMetaFile = fullfile(rootDir, 'globalBorder.mat');
p.svg.graphFile = fullfile(rootDir, 'graph.mat');
p.svg.correspondenceFile = fullfile(rootDir, 'correspondences.mat');
p.svg.aggloPredFile = fullfile(rootDir, 'segmentAggloPredictions.mat');
p.svg.heuristicFile = fullfile(rootDir, 'heuristicResultComplete.mat'); % assumes scores exist for each segment but not available yet

p.connectome = struct;
p.connectome.saveFolder = fullfile(rootDir, 'connectome');

p.agglo = struct;
p.agglo.axonAggloFile = axonAggloFile;
p.agglo.spineHeadFile = spineHeadFile;

axons = load(p.agglo.axonAggloFile, 'axons');
axons = axons.axons;

axons = arrayfun( ...
    @(a) double(a.segIds), axons, ...
    'UniformOutput', false);

%% interface agglomeration parameters
params_agglo = struct;
params_agglo.synT = [];
params.params_agglo = params_agglo;

% NOTE(amotta): See
% L23.Synapse.Script.plotPrecisionRecallScores
% git@gitlab.mpcdf.mpg.de:connectomics/amotta.git a5b7cdd5165ad024ad08d9705987cef6c17cb590
% amotta@m-01522. MATLAB 9.3.0.713579 (R2017b). 12-Dec-2019 19:08:37
params.synT_spine = -1.6;

% NOTE(amotta): Let's use the same thresholds for shafts, for now. We don't
% have test set for shaft synapes, so...
params.synT_shaft = params.synT_spine;
                             
%% determine interfaces onto spines
m = load(p.svg.edgeFile, 'edges');
edges = m.edges;
maxEdgeIdx = size(edges, 1);
m = load(p.agglo.spineHeadFile, 'shAgglos');
isSHInt = any(ismember(edges, cell2mat(m.shAgglos)), 2);

%% use separate scores for spine/other interfaces
scores_sh = SynEM.Util.loadGlobalSynScores(p.svg.synScoreFile, maxEdgeIdx);
scores = scores_sh;

% determine initial synaptic indices (separately for spine/other)
synIdx = false(size(scores, 1), 1);
synIdx(isSHInt) = any(scores(isSHInt, :) > params.synT_spine, 2);
synIdx(~isSHInt) = any(scores(~isSHInt, :) > params.synT_shaft, 2);

% convert scores to graph (i.e. including correspondences)
m = load(p.svg.graphFile, 'borderIdx');
borderIdx = m.borderIdx;
tmp = nan(size(borderIdx, 1), 2);
tmp(~isnan(borderIdx), :) = scores;
scores = tmp;

% convert synIdx to graph
tmp = false(size(borderIdx, 1), 1);
tmp(~isnan(borderIdx)) = synIdx;
synIdx = tmp;

%% interface agglomeration
[synapses, isSpineSyn, heuristics, debug] = ...
    L4.Synapses.createSynapseAgglomerates( ...
        p, params_agglo, axons, [], scores, synIdx);

%{
%% get SVM predictions to correct synapse agglomerates
boxIds = 1:prod(p.tiles);
tic;
endStep = numel(boxIds);
for boxId = 1:endStep
    % load SVM predictions
    m = load([p.local(boxId).saveFolder 'svmSegData.mat'], seg, synCom, vcCom, miCom);

end  
%}
Util.save( ...
    outputFile, synapses, isSpineSyn, ...
    heuristics, debug, params, info);
