function job = restartTasksWithErrors(job, cluster)
% Pass job object, will restart all tasks that had errors as new job, pass 2nd argument is cluster object to be used
    
    errorCell = {job.Tasks(:).Error};
    idxError = ~cellfun(@isempty, errorCell);
    inputCell = {job.Tasks(idxError).InputArguments};
    % This is now necessary with Cluster.startJob, any way to solve better?
    inputCell = cellfun(@(x){x}, inputCell, 'uni', 0);
    functionH = job.Tasks(find(idxError,1)).Function;
    job = Cluster.startJob(functionH, inputCell, 'name', 'restartedTasks', 'cluster', cluster);
end

