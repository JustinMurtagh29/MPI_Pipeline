function job = startCPU(fH, iC, jN, requiredMemory, group, priority, rt)
    % Wrapper function for startJob.m used for backward compability
    % Set default values for additional input arguments
    if ~exists('requiredMemory', 'var') || isempty(requiredMemory)
        requiredMemory = 12;
    end
    if ~exists('group', 'var') || isempty(group)
        group = 1;
    end
    if ~exists('priority', 'var') || isempty(priority)
        priority = -500;
    end
    if ~exists('rt', 'var') || isempty(rt)
        rt = 29;
    end
    
    % Translate parameters for `Cluster` module
    memory = ceil(requiredMemory);
    priority = ceil((priority + 1000) / 10);
    time = sprintf('%02d:%02d:00', floor(rt), ceil(60 * mod(rt, 1)));

    job = Cluster.startJob( ...
        fH, iC, ...
        'name', jN, ...
        'cluster', { ...
            'time', time, ...
            'memory', memory, ...
            'priority', priority}, ...
        'taskGroupSize', group);
end

