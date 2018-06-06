% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
%   Benedikt Staffler <benedikt.staffler@brain.mpg.de>
%
% Differences to connectEM.Axon.Script.calculateSynToSynDists:
%   - new synapse agglo version and the corresponding connectome

clear;


%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_spine_attachment.mat');
connFile = fullfile(rootDir, 'connectomeState', ...
    'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v6-somaH.mat');
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos-v6-somaH.mat');

outDir = fullfile(rootDir, 'connectomeState');

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

if ~exist(outFile, 'file')
    Util.saveStruct(outFile, out);
    Util.protect(outFile);
    Util.log('Saving output to %s.', outFile);
else
    Util.log('File %s already exists.', outFile);
end
