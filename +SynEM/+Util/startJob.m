function job = startJob( cluster, fh, inputCell, fhArgsOut, jobName )
%STARTJOB Job creation wrapper.
% INPUT cluster: A parcluster object.
%       fh: Handle to the function that is called by the job.
%       inputCell: Cell array of cell arrays. The outer cell array
%           corresponds to the number of tasks of the job. Each outer cell
%           contains a cell array with all inputs to fh for one task.
%       fhArgsOut: (Optional) Number of output arguments from fh.
%                  (Default: 0)
%       jobName: (Optional) Name for the job.
%                (Default: Name set by MATLAB).
% OUTPUT job: A job object.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%create job object
job = createJobPwd(cluster);
if exist('jobName','var') && ~isempty(jobName)
    job.Name = jobName;
end

if ~exist('fhArgsOut','var') || isempty(fhArgsOut)
    fhArgsOut = 0;
end

createTask(job, fh, fhArgsOut, inputCell, 'CaptureDiary', true);

submit(job);

end

