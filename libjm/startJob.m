function job = startJob(jm, functionHandle, inputCell, jobName, group )

    if nargin < 5
        group = 1;
    end
    % Create job on cluster
    job = createJob(jm);
    % Pass path for now instead of attatching files, otherwise passing anonymous functions will not work!!!
    ap = strsplit(path, ':');
    job.AdditionalPaths = ap;
    job.AutoAttachFiles = false;
    % Name job
    job.Name = [jobName '_' datestr(clock, 30)];
    if group ~= 1
        % Submit tasks together as batches (use in case single task is too short to avoid SGE/matlab starting overhead)
        grouping = 0:group:length(inputCell);
        if grouping(end) ~= length(inputCell)
            grouping(end+1) = length(inputCell);
        end
        % Group together tasks
        fH = @jobWrapper;
        iC = cell(length(grouping)-1,1);
        for i=1:length(grouping)-1
            iC{i} = {functionHandle, inputCell(grouping(i)+1:grouping(i+1))};
        end
        % Create regrouped task(s)
        createTask(job, fH, 0, iC, 'CaptureDiary', true);
    else
        % Create Task(s)
        createTask(job, functionHandle, 0, inputCell, 'CaptureDiary', true);
    end
    % Start job
    submit(job);

end

