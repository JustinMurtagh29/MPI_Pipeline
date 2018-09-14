%{ 
author: Yagmur Yener
email: yagmur.yener.yy@gmail.com

run locally

Locating pre- and post- synapses in the astrocyte annotated segmented
volume
%}

syn = load('~/GABA/astrocyte/synapses/syn.mat');
syn = syn.syn;
seg = load('~/GABA/astrocyte/synapses/seg.mat');
seg = seg.seg;
conn = load('~/GABA/astrocyte/synapses/conn.mat');
conn = conn.conn;
shAgglos = load('~/GABA/astrocyte/synapses/shAgglos.mat');
shAgglos = shAgglos.shAgglos;
borderAreas = load('~/GABA/astrocyte/synapses/borderAreas.mat');
borderAreas = borderAreas.borderAreas;