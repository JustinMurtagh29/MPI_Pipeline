function myelinCount = collectMyelinData(p)


inputCell = num2cell(1:numel(p.local));
cluster = Cluster.getCluster( ...
    '-pe openmp 1', ...
    '-p -0', ...
    '-l h_vmem=12G', ...
    '-l s_rt=4:00:00', ...
    '-l h_rt=4:00:30');
job = Cluster.startJob(@Myelin.collectMyelinDataCube, inputCell,'sharedInputs',{p}, 'cluster', cluster, 'name', 'myelinDataCollection', 'taskGroupSize',1 ,'numOutputs',1);
Cluster.waitForJob(job);

outp = Cluster.fetchTaskOutputs(job);


myelinCount = cat(1,outp{:});

Util.save(fullfile(p.saveFolder,'globalMyelinCount.mat'),myelinCount)