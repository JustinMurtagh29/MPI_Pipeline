%{ 
author: Yagmur Yener
email: yagmur.yener.yy@gmail.com

run locally

unique map from synapses to axons
%}

% syn = load('~/GABA/astrocyte/synapses/syn.mat');
% syn = syn.syn;
% seg = load('~/GABA/astrocyte/synapses/seg.mat');
% seg = seg.seg;
conn = load('~/GABA/astrocyte/synapses/conn.mat');
conn = conn.conn;


%%

syn_id_reps = cellfun(@(id_arrays) numel(id_arrays), conn.connectome.synIdx);
axon_ids_rptd = repelem(conn.connectome.edges(:,1), syn_id_reps);

id_syn_axon = [cell2mat(conn.connectome.synIdx) axon_ids_rptd];

%%

lut_syn_axon = zeros(max(max(id_syn_axon(:,1))),1);
lut_syn_axon(id_syn_axon(:,1)) = id_syn_axon(:,2);

%%
syn_axon_segs = cell(size(lut_syn_axon));

syn_axon_segs(id_syn_axon(:,1)) = conn.axons(id_syn_axon(:,2));

%% 
save('~/GABA/astrocyte/synapses/synID_axonSeg.mat', 'syn_axon_segs')