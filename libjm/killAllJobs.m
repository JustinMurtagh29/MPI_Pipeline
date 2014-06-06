function killAllJobs(which,username, finished)
%killJob( jobIDsToKill )
%   Kills jobs on "which" jobmanager (if finsished = 1, only finshed jobs)

jm = findJm;
if any(which == 1)
    if finished
        jobs = findJob(jm(1), 'Username', username, 'state', 'finished');
    else
        jobs = findJob(jm(1), 'Username', username);
    end
    if ~isempty(jobs)
        destroy(jobs);
    end
end
if any(which == 2)
    if finished
        jobs = findJob(jm(2), 'Username', username, 'state', 'finished');
    else
        jobs = findJob(jm(2), 'Username', username);
    end
    if ~isempty(jobs)
        destroy(jobs);
    end 
end


end

