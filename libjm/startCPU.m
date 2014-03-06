function job = startCPU( functionHandle, inputCell )

% Load cluster configuration
jm = findResource('scheduler', 'type', 'jobmanager', 'configuration', 'fermat-cpu');
% Create job on cluster
job = createJob(jm, 'configuration', 'fermat-cpu');
createTask(job, functionHandle, 0, inputCell, 'configuration', 'fermat-cpu');
% Start job
submit(job);

end

