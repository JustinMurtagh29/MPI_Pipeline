function job = startCPU( functionHandle, inputCell, jobName )
    global GLOBAL_HOST;
    global GLOBAL_CODE_DIR;

    pathDependencies = genpathGit(GLOBAL_CODE_DIR);
    pathDependencies = regexp(pathDependencies, ':', 'split');
    pathDependencies = pathDependencies(cellfun(@(x)not(isempty(x)), pathDependencies));

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
        job.AdditionalPaths = genpathGit(GLOBAL_CODE_DIR); 
        for i=1:length(functionHandle)
            createTask(job, functionHandle{i}, 0, inputCell{i}, 'MaximumNumberOfRetries', 5, 'Timeout', 90000, 'CaptureCommandWindowOutput', 0);
        end
    end

    % Start job
    submit(job);

end

