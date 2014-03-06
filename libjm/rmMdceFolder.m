function rmMdceFolder(jm, dir)
%% Get ownership of folder owned by mdceproc
job = createJob(jm);
task = createTask(job, @system, 0, {['rm -rf ' dir]});
submit(job);
end

