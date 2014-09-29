function job = startCPU( functionHandle, inputCell, jobName )
global GLOBAL_HOST;

pathDependencies = {'/home/trende/august/' '/home/trende/august/libjm'};

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

