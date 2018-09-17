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
graph = load('~/GABA/astrocyte/synapses/graph.mat');
graph = graph.graph;
maxSegId = 15030572; %maximum possible segment ID

%%
asiT = connectEM.Connectome.buildAxonSpineInterfaces( ...
        maxSegId, graph, shAgglos, conn, syn);
%%

save('~/GABA/astrocyte/synapses/asiT.mat', 'asiT')