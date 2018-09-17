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
% Segment ids to Presynapses mapping
lut_seg_presyn = Agglo.buildLUT(maxSegId, syn.synapses.presynId);

numel(intersect(find(lut_seg_sh) , find(lut_seg_presyn) )) %segments shared by synapses and spine heads
numel(intersect(intersect(find(lut_seg_sh) , find(lut_seg_presyn) ), find(lut_seg))) %segments shared by synapses and spine heads in this vol
%%
lut_presyn_sh = cell(size(syn.synapses.presynId));
syn_idx = find(lut_syn);
sh_idx = find(lut_sh);
for i = syn_idx
    preSynArray = syn.synapses.presynId{i};
    
    temp = cellfun(@(x) any(x==preSynArray),  shAgglos);
    lut_presyn_sh{i} = find(temp);
    
    
end

preSynSegments = syn.synapses.presynId;
%tell if the synapse id in the box given its total segments
lut_sh_idx = cellfun(@(preSegIds) any(shAgglos==preSegIds), preSynSegments); %look-up table for syns
