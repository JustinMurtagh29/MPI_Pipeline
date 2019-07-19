%{ 
author: Yagmur Yener
email: yagmur.yener.yy@gmail.com

run on gaba server

Saving the synapse information to files for faster usage locally.
%}

%% 

rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);
conn = load(connFile); %axons, dendrites, connectome
syn = load(conn.info.param.synFile); % synapses (segment IDs)


%%
save('/gaba/u/yyener/astrocyte/synapses/partiallySplit/syn.mat', 'syn')
save('/gaba/u/yyener/astrocyte/synapses/partiallySplit/param.mat', 'param')
save('/gaba/u/yyener/astrocyte/synapses/partiallySplit/conn.mat', 'conn')




