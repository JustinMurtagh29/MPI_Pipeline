function jobID = extractJobId(cmdOut)
% Extracts the job ID from the qsub command output for SGE

% Copyright 2010-2011 The MathWorks, Inc.

% The output of qsub will be:
% Your job 496 ("<job name>") has been submitted

% Now parse the output of bsub to extract the job number
jobNumberStr = regexp(cmdOut, 'job [0-9]*', 'once', 'match');
jobID = sscanf(jobNumberStr, 'job %d');
dctSchedulerMessage(0, '%s: Job ID %d was extracted from qsub output %s.', mfilename, jobID, cmdOut);

