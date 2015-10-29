function job = startJob(jm, functionHandle, inputCell, jobName )
    global GLOBAL_CODE_DIR;

    % Get and format path dependencies  
    pathDependencies = genpathGit(GLOBAL_CODE_DIR);
    pathDependencies = regexp(pathDependencies, ':', 'split');
    pathDependencies = pathDependencies(cellfun(@(x)not(isempty(x)), pathDependencies));
    % Create job on cluster
    job = createJob(jm);
    job.AutoAttachFiles = false;
    job.Name = [datestr(clock, 30) '_' jobName];
    job.AdditionalPaths = pathDependencies; 
    task = createTask(job, functionHandle, 0, inputCell, 'CaptureDiary', true);
    % Start job
    submit(job);

end

