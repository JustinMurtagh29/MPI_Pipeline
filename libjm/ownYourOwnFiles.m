function ownYourOwnFiles(jm, directory)
%% Get ownership (was only possible as root) or at least r+w access (mdceproc) to a folder
job = createJob(jm);
task = createTask(job, @system, 0, {['chmod -R 770 ' directory]});
submit(job);
end

