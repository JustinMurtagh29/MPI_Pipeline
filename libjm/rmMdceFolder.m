function rmMdceFolder(dir)
% Remove folder owned by mdceproc
jm = findJm;
jm = jm(1);
job = createJob(jm);
task = createTask(job, @system, 0, {['rm -rf ' dir]});
submit(job);
end

