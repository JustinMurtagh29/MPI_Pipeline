function job = startGPU( functionHandle, inputCell, jobName )

% We only need CNN path dependencies on the GPU part of the cluster
pathDependencies = {'/usr/local/jacket/' '/usr/local/jacket/engine/' '/zdata/manuel/code/CNN/' '/zdata/manuel/code/auxiliaryMethods/' '/zdata/manuel/code/auxiliaryMethods/cubes/'};

% Load cluster configuration
jm = findJm();
% Create job on cluster
% 'PathDependencies', pathDependencies, 
job = createJob(jm(2), 'RestartWorker', 1, 'PathDependencies', pathDependencies, 'Name', jobName);
for i=1:length(functionHandle);
	createTask(job, functionHandle{i}, 0, inputCell{i}, 'MaximumNumberOfRetries', 5, 'Timeout', 90000, 'CaptureCommandWindowOutput', 0);
end
% Start job
submit(job);

end

