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
maxSegId = 15030572; %maximum possible segment ID

graph = 1;

%%
asiT = buildAxonSpineInterfaces( ...
        param, graph, shAgglos, conn, syn)
