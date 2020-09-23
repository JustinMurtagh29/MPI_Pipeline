% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
connFile = fullfile(rootDir, 'connectome', 'Connectome_20191227T220548-results_20191227T220548-results-auto-spines-v3_SynapseAgglomerates--20191227T220548-results--20191227T220548-results-auto-spines-v3--v1.mat');
synScoreFile = fullfile(rootDir, 'globalSynapseScores.mat');

outVersion = 1;

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = load(connFile);

synFile = conn.info.param.synFile;
syn = load(synFile);

shAgglos = syn.info.param.spineHeadFile;
shAgglos = Util.load(shAgglos, 'shAgglos');

somaAgglos = conn.denMeta.targetClass == 'Somata';
somaAgglos = conn.dendrites(somaAgglos);

%% Load synapse scores per edge
% This is based on Seg.IO.loadGraph
curEdgeToBorderIds = Util.load( ...
    fullfile(rootDir, 'graph.mat'), 'borderIdx');
curEdgeCount = numel(curEdgeToBorderIds);
curEdgeToBorderIds = find(not(isnan(curEdgeToBorderIds)));

[curSynEdgeIds, curSynScores] = ...
    Util.load(synScoreFile, 'edgeIdx', 'synScores');
curSynEdgeIds = curEdgeToBorderIds(curSynEdgeIds);

synScores = nan(curEdgeCount, 2);
synScores(curSynEdgeIds, :) = curSynScores;

%% Compute synapse types
clear cur*;

synTypes = connectEM.Synapse.classifyType( ...
    param, syn.synapses, synScores, shAgglos, somaAgglos);

%% Save result
out = struct;
out.info = info;
out.types = synTypes;

[outDir, outFile] = fileparts(synFile);
outFile = sprintf('%s__types_v%d.mat', outFile, outVersion);
outFile = fullfile(outDir, outFile);

Util.saveStruct(outFile, out);
Util.protect(outFile);
