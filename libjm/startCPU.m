function job = startCPU(fH, iC, jN, requiredMemory, group, priority)
    % Wrapper function for startJob.m used for backward compability

    % Set default values for additional input arguments
    if nargin < 4
        requiredMemory = 12;
    end
    if nargin < 5
        group = 1;
    end
    if nargin < 6
        priority = -200;
    end

    clusterCPU = Cluster.getCluster( ...
        '-pe openmp 1', ...
        ['-p ' num2str(priority)], ...
        ['-l h_vmem=' num2str(requiredMemory) 'G'], ...
        '-l s_rt=23:50:00', ...
        '-l h_rt=24:00:00');

    job = Cluster.startJob(fH, iC, 'cluster', clusterCPU, 'name', jN, 'taskGroupSize', group);

end

