function job = startCPU( functionHandle, inputCell, jobName )
global GLOBAL_HOST;

pathDependencies = {'/zdata/manuel/code/CNN/' '/zdata/manuel/code/auxiliaryMethods/' '/zdata/manuel/code/auxiliaryMethods/cubes/', ...
	'/zdata/manuel/code/active/' '/zdata/manuel/code/correspondence/' '/zdata/manuel/code/graphConstruction/', ...
	'/zdata/manuel/code/filterbank/'};

% Load cluster configuration
jm = findJm();
% Create job on cluster
if strcmp(GLOBAL_HOST,'fermat01')
    job = createJob(jm(1), 'RestartWorker', true, 'PathDependencies', pathDependencies, 'Name', jobName);
    for i=1:length(functionHandle)
	    createTask(job, functionHandle{i}, 0, inputCell{i}, 'MaximumNumberOfRetries', 5, 'Timeout', 90000, 'CaptureCommandWindowOutput', 0);
    end
elseif strcmp(GLOBAL_HOST, 'gaba')
    job = createJob(jm(1), 'RestartWorker', true, 'Name', jobName);
    job.AdditionalPaths = {'/gaba/u/mberning/code/CNN/'}; 
    for i=1:length(functionHandle)
	    createTask(job, functionHandle{i}, 0, inputCell{i});
    end
end

% Start job
submit(job);

end

