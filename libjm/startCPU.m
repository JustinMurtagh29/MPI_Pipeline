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
        priority = -500;
    end
    
    cluster = parcluster('local');
    job = Cluster.startJob(fH, iC, 'cluster', clusterCPU, 'name', jN, 'taskGroupSize', group);

end

