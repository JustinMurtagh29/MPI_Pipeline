function tests = testJobSubmission()
    tests = functiontests(localfunctions);
end

function testCaseSimple(testResult)
    % Test simple job with one output argument
    functionHandle = @times;
    inputCell = arrayfun(@(x){{randn(x), randn(x)}}, 1:10);
    
    % Run once on local cluster and once directly
    resultJob = runJobAndReturnResult(functionHandle, inputCell);
    resultLocal = cellfun(@(x)functionHandle(x{:}), inputCell, 'uni', 0);
    
    % Compare results (more likely error thrown above if sth. goes wrong)
    % cellfun and job submission yield tranposed cell arrays, fix?
    verifyEqual(testResult, resultJob, resultLocal');
end

function testCaseSharedInputs(testResult)
    % Test simple job with one output argument and shared inputs
    functionHandle = @(x,y,z)times(x,y).^z;
    inputCell = arrayfun(@(x){{randn(x)}}, 1:10);
    sharedInputs = {randn(1) rand(1)};
    sharedInputsLocation = [1 3];
    inputs = cellfun(@(x) [sharedInputs(1), x(:), sharedInputs(2)], ...
        inputCell, 'uni', 0);
    
    % Run once on local cluster and once directly
    resultJob = runJobAndReturnResult(functionHandle, inputCell, ...
        'sharedInputs', sharedInputs, ...
        'sharedInputsLocation', sharedInputsLocation);
    resultLocal = cellfun(@(x)functionHandle(x{:}), inputs, 'uni', 0);
    
    % Compare results (more likely error thrown above)
    verifyEqual(testResult, resultJob, resultLocal');
end

function testCaseGrouping(testResult)
    % Test simple job with one output argument and grouped tasks
    functionHandle = @(x,y,z)times(x,y).^z;
    inputCell = arrayfun(@(x){{randn(1) randn(x) randn(1)}}, 1:10);
    taskGroupSize = 3;
    
    % Run once on local cluster and once directly
    resultJob = runJobAndReturnResult(functionHandle, inputCell, ...
        'taskGroupSize', taskGroupSize);
    resultLocal = cellfun(@(x)functionHandle(x{:}), inputCell, 'uni', 0);
    
    % See Cluster.startJob documentation, outputs needs reshaping for jobs
    % that use taskGroupSize
    resultJob = vertcat(resultJob{:});
    % Compare results (more likely error thrown above)
    verifyEqual(testResult, resultJob, resultLocal');
end

function resultJob = runJobAndReturnResult(functionHandle, inputCell, varargin)
    % Run function using generic local cluster
    cluster = parcluster();
    job = Cluster.startJob(functionHandle, inputCell, ...
        'cluster', cluster, 'numOutputs', 1, varargin{:});
    wait(job, 'finished');
    resultJob = fetchOutputs(job);
    delete(job);
end
