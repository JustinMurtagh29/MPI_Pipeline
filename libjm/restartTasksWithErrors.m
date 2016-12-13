function job = restartTasksWithErrors(job, gpuFlag)
% Pass job object, will restart all tasks that had errors as new job, pass 2nd argument true if GPU is required
    
    if nargin == 1
        gpuFlag = false;
    end
    errorCell = {job.Tasks(:).Error};
    idxError = ~cellfun(@isempty, errorCell);
    inputCell = {job.Tasks(idxError).InputArguments};
    functionH = job.Tasks(find(idxError,1)).Function;
    if gpuFlag
        job = startGPU(functionH, inputCell, 'restartedTasks');
    else
        job = startCPU(functionH, inputCell, 'restartedTasks');
    end
end
