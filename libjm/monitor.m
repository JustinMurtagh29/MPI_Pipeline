function monitor()

jmGPU = findResource('scheduler', 'type', 'jobmanager', 'configuration', 'fermat-cnn');
jmCPU = findResource('scheduler', 'type', 'jobmanager', 'configuration', 'fermat-cpu');

while true
	clc;
	display('Jobmanager for CPU:');
	display(['Number Workers:' 9 9  num2str(jmCPU.ClusterSize, '%.3i')]);
	display(['Busy Workers:' 9 9 num2str(jmCPU.NumberOfBusyWorker, '%.3i')]);
	display(['Idle Workers:' 9 9 num2str(jmCPU.NumberOfIdleWorker, '%.3i')]);
	display('Jobmanager for GPU:');
	display(['Number Workers:' 9 9  num2str(jmGPU.ClusterSize, '%.3i')]);
	display(['Busy Workers:' 9 9 num2str(jmGPU.NumberOfBusyWorker, '%.3i')]);
	display(['Idle Workers:' 9 9 num2str(jmGPU.NumberOfIdleWorker, '%.3i')]);
	gpuJobs = findJob(jmGPU, 'state', 'running');
	cpuJobs = findJob(jmCPU, 'state', 'running');
	display(cpuJobs);
	display(gpuJobs);
	pause(5*60);
end

end
