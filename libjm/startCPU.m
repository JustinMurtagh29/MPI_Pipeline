function job = startCPU(fH, iC, jN, requiredMemory, group, priority, rt);
    % Wrapper function for startJob.m used for backward compability
    if nargin < 7 || isempty(rt)
	rt = 23;
    end
    % Set default values for additional input arguments
    if nargin < 4 || isempty(requiredMemory)
        requiredMemory = 12;
    end
    if nargin < 5 || isempty(group)
        group = 1;
    end
    if nargin < 6 || isempty(priority)
        priority = -500;
    end

    clusterCPU = Cluster.getCluster( ...
        '-pe openmp 1', ...
        ['-p ' num2str(priority)], ...
        ['-l h_vmem=' num2str(requiredMemory) 'G'], ...
	sprintf('-l s_rt=%02d:%02d:00',floor(rt),round((rt-floor(rt))*60)), ...
	sprintf('-l h_rt=%02d:%02d:30',floor(rt),round((rt-floor(rt))*60)));

    job = Cluster.startJob(fH, iC, 'cluster', clusterCPU, 'name', jN, 'taskGroupSize', group);

end

