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

%% Load the validation results from file

%location in the L4 dataset
box_offset = [4179, 4994, 2264];
box_shape = [178, 178, 72];
bbox = [box_offset' , (box_offset+box_shape-1)'];
%%
%Region annotated by CNN (voxel-wise)
%pred: binary thresholded predictions 1 means astro
%prediction: predicted probabilities continuous between -1.7 and 1.7
%(scaled tanh output). -1 means astro
astro_annot = load('/gaba/u/yyener/astrocyte/predictions/unet_aug/v4_val.mat');

%% Get the volumes for synapses and astrocytes

syn_points = Seg.Global.getSegToPointMap(param); %segments to point indices


%% Save to a mat file for local processing
save('/gaba/u/yyener/astrocyte/synapses/syn_points.mat', 'syn_points')

%%
save('/gaba/u/yyener/astrocyte/synapses/syn.mat', 'syn')
save('/gaba/u/yyener/astrocyte/synapses/param.mat', 'param')
save('/gaba/u/yyener/astrocyte/synapses/conn.mat', 'conn')

%% segmentation volumes

% for speed
param.seg = struct;
param.seg.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
param.seg.backend = 'wkwrap';
%save('/gaba/u/yyener/astrocyte/synapses/param.mat', 'param')

seg = loadSegDataGlobal(param.seg, bbox);
save('/gaba/u/yyener/astrocyte/synapses/seg.mat', 'seg')



