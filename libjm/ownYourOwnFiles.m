function ownYourOwnFiles(directory)
%% Get ownership (was only possible as root) or at least r+w access (mdceproc) to a folder
jm = findJm();
job = createJob(jm(1));
task = createTask(job, @system, 0, {['chmod -R 770 ' directory]});
submit(job);
end

