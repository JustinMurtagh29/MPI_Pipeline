function isReadable = testFolderAvailibility( path )
jm = findJm();
job = createJob(jm(1));
task = createTask(job, @exist, 1, {path 'file'});
submit(job);
waitForState(job, 'finished');
isReadable = getAllOutputArguments(job);
end

