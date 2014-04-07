function isReadable = testFolderAvailibility( path )
jm = findJm();
job = createJob(jm(2));
for i=1:28
    createTask(job, @exist, 1, {path 'file'});
end
submit(job);
waitForState(job, 'finished');
isReadable = getAllOutputArguments(job);
end

