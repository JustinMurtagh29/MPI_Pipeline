function restartTasksWithErrors(job, gpuFlag)
% Pass job object, will restart all tasks that had errors as new job, pass 2nd argument true if GPU is required
    
    if nargin == 1
        gpuFlag = false;
    end
    errorCell = {job.Tasks(:).Error};
    idxError = cellfun(@isempty, errorCell);
    inputCell = {job.Tasks(~idxError).InputArguments};
    functionH = @galleryCortex;
    if gpuFlag
        job = startGPU(functionH, inputCell, 'restarted tasks');
    else
        job = startCPU(functionH, inputCell, 'restarted tasks');
    end
end

