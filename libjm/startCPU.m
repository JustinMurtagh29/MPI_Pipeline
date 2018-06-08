function job = startCPU(fH, iC, jN, requiredMemory, group, priority, rt);
    % Wrapper function for startJob.m used for backward compability
    if nargin < 7 || isempty(rt)
	rt = 29;
    end
    % Set default values for additional input arguments
    if nargin < 4 || isempty(requiredMemory)
        requiredMemory = 12;
    end
    if nargin < 5 || isempty(group)
        group = 1;
    end
    if nargin < 6 || isempty(priority)
        priority = 50;
    end
    
    clusterCPU = Cluster.config( ...
        'priority',priority,...
        'memory', requiredMemory, ...
	    'time',[num2str(floor(rt),'%02d') ':' num2str(round((rt-floor(rt))*60),'%02d') ':' num2str(0,'%02d')]);
    job = Cluster.startJob(fH, iC, 'cluster', clusterCPU, 'name', jN, 'taskGroupSize', group);

end

