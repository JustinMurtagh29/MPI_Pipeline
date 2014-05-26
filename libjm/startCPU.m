function startCPU( functionHandle, inputCell, jobName )

pathDependencies = {'/zdata/manuel/code/CNN/' '/zdata/manuel/code/auxiliaryMethods/' '/zdata/manuel/code/auxiliaryMethods/cubes/', ...
       '/zdata/manuel/code/active/' '/zdata/manuel/code/correspondence/' '/zdata/manuel/code/graphConstruction/'};

% Load cluster configuration
jm = findJm();
jm = jm(1);
% Create job on cluster
job = createJob(jm, 'RestartWorker', 1, 'PathDependencies', pathDependencies, 'Name', jobName);
for i=1:length(functionHandle)
	createTask(job, functionHandle{i}, 0, inputCell{i}, 'MaximumNumberOfRetries', 5, 'Timeout', 90000, 'CaptureCommandWindowOutput', 0);
end
% Start job
submit(job);

end

