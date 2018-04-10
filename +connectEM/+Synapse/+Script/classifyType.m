% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_spine_attachment.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');

outDir = '/home/amotta/Desktop';

[~, outFile] = fileparts(synFile);
outFile = sprintf('%s_classified.mat', outFile);
outFile = fullfile(outDir, outFile);

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

%% Generate output
out = syn;
out.info = info;

Util.saveStruct(outFile, out);
Util.protect(outFile);
