%{ 
author: Yagmur Yener
email: yagmur.yener.yy@gmail.com

run locally

Completing the synapses using segments from spine heads
%}

syn = load('~/GABA/astrocyte/synapses/syn.mat');
syn = syn.syn;
seg = load('~/GABA/astrocyte/synapses/seg.mat');
seg = seg.seg;
conn = load('~/GABA/astrocyte/synapses/conn.mat');
conn = conn.conn;
shAgglos = load('~/GABA/astrocyte/synapses/shAgglos.mat');
shAgglos = shAgglos.shAgglos;
%%

maxSegId = 15030572; %maximum possible segment ID

%Segments in this Volume
lut_seg = false(maxSegId, 1); %initialize as false
lut_seg(setdiff(seg, 0)) = true; %sets sed ids true if unique seg nonzero

% Synapses in this volume
%merge segments of pre and post
synSegments = cellfun(@vertcat, syn.synapses.presynId, syn.synapses.postsynId, 'UniformOutput', false);
%tell if the synapse id in the box given its total segments
lut_syn = cellfun(@(synSegIds) all(lut_seg(synSegIds) == true), synSegments); %look-up table for syns

% Spine heads in this volume
lut_sh = cellfun(@(shSegId) all(lut_seg(shSegId) == true), shAgglos);

%%

% Segment ids to Spine Heads mapping
lut_seg_sh = Agglo.buildLUT(maxSegId, shAgglos);
% Segment ids to Postsynapses mapping
lut_seg_postsyn = Agglo.buildLUT(maxSegId, syn.synapses.postsynId);
% Segment ids to Presynapses mapping
lut_seg_presyn = Agglo.buildLUT(maxSegId, syn.synapses.presynId);
% Segment ids to axons mapping
lut_seg_axon = Agglo.buildLUT(maxSegId, conn.axons);

numel(intersect(find(lut_seg_sh) , find(lut_seg_presyn) )) %segments shared by Presynapses and spine heads (647)
numel(intersect(find(lut_seg_sh) , find(lut_seg_postsyn) )) %segments shared by Postsynapses and spine heads (268601)
numel(intersect(intersect(find(lut_seg_sh) , find(lut_seg_postsyn) ), find(lut_seg))) %segments shared by synapses and spine heads in this vol

%% Complete presynaptic side with spine heads

common_segs = intersect(intersect(find(lut_seg_sh) , find(lut_seg_postsyn) ), find(lut_seg));
postSynCompleted = syn.synapses.postsynId;
for i = 1:numel(common_segs)

sh_id = lut_seg_sh(common_segs(i));
syn_id = lut_seg_postsyn(common_segs(i));

postSynCompleted{syn_id} = shAgglos{sh_id};
end

%% Complete postsynaptic side with axons








