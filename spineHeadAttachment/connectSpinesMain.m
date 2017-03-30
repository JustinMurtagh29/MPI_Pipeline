function connectSpinesMain(makeConnM)
if makeConnM
    graph = load('/gaba/u/mberning/results/pipeline/20151111T183414/graph.mat', 'edges', 'prob');
    result = agglomerate_fast_makeM(graph, 0, {}, [], graph.prob);
    connM = result.connM + result.connM';
    save('/gaba/u/kboerg/connM.mat', 'connM');
else
addpath('/gaba/u/kboerg/auxiliaryMethods');
cluster = Cluster.getCluster('-pe openmp 1 -l h_vmem=6G,h_rt=100:00:00 -p -500');
job = createJob(cluster);
for i = 1 : 878
   inputCell{i} = {i};
end
createTask(job, @connectSpines, 0, inputCell)
submit(job);
end