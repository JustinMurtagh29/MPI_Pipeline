% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

thisDir = fileparts(mfilename('fullpath'));
ctrlDir = fullfile(thisDir, 'annotations', 'random-spine-synapses');

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

%% Loading augmented graph
graph = Graph.load(rootDir);
graph(~graph.borderIdx, :) = [];

borderAreas = fullfile(rootDir, 'globalBorder.mat');
borderAreas = load(borderAreas, 'borderArea2');
borderAreas = borderAreas.borderArea2;

graph.borderArea = borderAreas(graph.borderIdx);
clear borderAreas;

%% Loading control data
ctrlFiles = dir(fullfile(ctrlDir, '*.nml'));
ctrlFiles = fullfile(ctrlDir, {ctrlFiles.name});

ctrlSynT = connectEM.Consistency.loadAnnotations(param, ctrlFiles);
ctrlSynT = vertcat(ctrlSynT{:});

% Calculate synapse areas
ctrlSynT.area = connectEM.Consistency.calcSynapseAreas(param, graph, ctrlSynT);
