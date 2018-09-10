%{ 
author: Yagmur Yener
email: yagmur.yener.yy@gmail.com

run on gaba server

The first test script for looking at the astrocyte annotated regions and
the synapse regions.
Trying on a small annotated region.
Saving the synapse information to files for faster usage locally.
%}

%% 

rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);
conn = load(connFile); %axons, dendrites, connectome
syn = load(conn.info.param.synFile); % synapses (segment IDs)

%% Get the volumes for synapses
% not the volume but indices of every point in a synapse segment

syn_points = Seg.Global.getSegToPointMap(param); %segments to point indices

%% Save to a mat file for local processing
save('/gaba/u/yyener/astrocyte/synapses/syn_points.mat', 'syn_points')






