function job = startCNN( cnet, stacks, settings )

% Load cluster configuration
jm = findResource('scheduler', 'type', 'jobmanager', 'configuration', 'fermat-cnn');

% Create Directories for saving results & start learning
if ~exist(cnet.run.savingPath, 'dir')
	mkdir(cnet.run.savingPath);
end

% Create job on cluster
job = createJob(jm, 'configuration', 'fermat-cnn');
inputargs = {cnet, stacks, settings};
createTask(job, @learn, 0, inputargs, 'configuration', 'fermat-cnn');
% Start job
submit(job);

end

