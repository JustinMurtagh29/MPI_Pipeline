function monitor(username)

jmGPU = findResource('scheduler', 'type', 'jobmanager', 'Username', username, 'configuration', 'fermat-cnn');
jmCPU = findResource('scheduler', 'type', 'jobmanager', 'Username', username, 'configuration', 'fermat-cpu');

while true
	clc;
	display('Jobmanager for CPU:');
	display(['Number Workers:' 9 9  num2str(jmCPU.ClusterSize, '%.3i')]);
	display(['Busy Workers:' 9 9 num2str(jmCPU.NumberOfBusyWorker, '%.3i')]);
	display(['Idle Workers:' 9 9 num2str(jmCPU.NumberOfIdleWorker, '%.3i')]);
    display(['-----------------------------------------------------------']);
    
    cpuJobs = findJob(jmCPU);    
    for i=1:length(cpuJobs)
        display(cpuJobs(i).Name);
        cpuTasks = findTask(cpuJobs(i), 'state', 'running');
	    display(['- Running tasks: '  9 num2str(length(cpuTasks))]);
        cpuTasks = findTask(cpuJobs(i), 'state', 'finished');
        display(['- Finished tasks: ' 9 num2str(length(cpuTasks))]);
        errorIdx = zeros(length(cpuTasks),1);
        for j=1:length(cpuTasks);
            errorIdx(j) = ~isempty(cpuTasks(j).Error);
        end
        display(['- Failed tasks: ' 9 num2str(sum(errorIdx))]);
    end 

    fprintf('\n\n\n\n');

    display('Jobmanager for GPU:');
	display(['Number Workers:' 9 9  num2str(jmGPU.ClusterSize, '%.3i')]);
	display(['Busy Workers:' 9 9 num2str(jmGPU.NumberOfBusyWorker, '%.3i')]);
	display(['Idle Workers:' 9 9 num2str(jmGPU.NumberOfIdleWorker, '%.3i')]);
	display(['-----------------------------------------------------------']);   

    gpuJobs = findJob(jmGPU);
    for i=1:length(gpuJobs)
	    display(gpuJobs(i).Name);
        gpuTasks = findTask(gpuJobs(i), 'state', 'running');
        display(['- Running tasks: ' 9 num2str(length(gpuTasks))]);
        gpuTasks = findTask(gpuJobs(i), 'state', 'finished');
        display(['- Finished tasks: ' 9 num2str(length(gpuTasks))]);
        errorIdx = zeros(length(gpuTasks),1);
        for j=1:length(gpuTasks);
            errorIdx(j) = ~isempty(gpuTasks(j).Error);
        end
        display(['- Failed tasks: ' 9 num2str(sum(errorIdx))]);

    end 
    pause(10);
end

end
