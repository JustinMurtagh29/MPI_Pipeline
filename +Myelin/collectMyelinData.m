function myelinCount = collectMyelinData(p)


inputCell = num2cell(1:numel(p.local));
cluster = Cluster.config( ...
    'priority', 0, ...
    'memory', 12, ...
    'time', '4:00:00');
job = Cluster.startJob(@Myelin.collectMyelinDataCube, inputCell,'sharedInputs',{p}, 'cluster', cluster, 'name', 'myelinDataCollection', 'taskGroupSize',1 ,'numOutputs',1);
Cluster.waitForJob(job);

outp = Cluster.fetchTaskOutputs(job);


myelinCount = cat(1,outp{:});

Util.save(fullfile(p.saveFolder,'globalMyelinCount.mat'),myelinCount)
