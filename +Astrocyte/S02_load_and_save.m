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

% maxSegId = Seg.Global.getMaxSegId(param);
% conn = load(connFile); %axons, dendrites, connectome
% syn = load(conn.info.param.synFile); % synapses (segment IDs)
% syn_points = Seg.Global.getSegToPointMap(param); %segments to point indices

%% Load the validation results from file

%location in the L4 dataset
box_offset = [4179, 4994, 2264];
box_shape = [178, 178, 72];
bbox = [box_offset' , (box_offset+box_shape-1)'];

% segmentation volumes

% for speed
param.seg = struct;
param.seg.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
param.seg.backend = 'wkwrap';

seg = loadSegDataGlobal(param.seg, bbox);
save('/gaba/u/yyener/astrocyte/synapses/seg.mat', 'seg')



%% Save to a mat file for local processing
% save('/gaba/u/yyener/astrocyte/synapses/syn_points.mat', 'syn_points')
% save('/gaba/u/yyener/astrocyte/synapses/syn.mat', 'syn')
% save('/gaba/u/yyener/astrocyte/synapses/param.mat', 'param')
% save('/gaba/u/yyener/astrocyte/synapses/conn.mat', 'conn')
% save('/gaba/u/yyener/astrocyte/synapses/param.mat', 'param')



