function job = startJob(jm, functionHandle, inputCell, jobName )

   % Create job on cluster
    job = createJob(jm);
    % Pass path for now instead of attatching files, otherwise passing anonymous functions will not work!!!
    ap = strsplit(path, ':');
    job.AdditionalPaths = ap;
    job.AutoAttachFiles = false;
    % Name job
    job.Name = [jobName '_' datestr(clock, 30)];
    % Create Task(s)
    task = createTask(job, functionHandle, 0, inputCell, 'CaptureDiary', true);
    % Start job
    submit(job);

end

