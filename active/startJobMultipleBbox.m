function job = startJobMultipleBbox( functionHandle, cnet, raw, class, bbox )

% Load cluster configuration
jm = findResource('scheduler', 'type', 'jobmanager', 'configuration', 'fermat-cpu');
% Create job on cluster
job = createJob(jm, 'configuration', 'fermat-cpu');
inputCell = {cnet, raw, class, bbox};
createTask(job, functionHandle, 0, inputCell, 'configuration', 'fermat-cpu');
% Start job
submit(job);

end

