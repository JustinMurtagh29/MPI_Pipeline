function job = createJobPwd(cluster)
%CREATEJOBPWD Create job and add current path with subdirectories to
%additional paths but disable AutoAttachFiles.
% INPUT cluster: A Matlab cluster object.
% OUTPUT job: A Matlab job object.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

p = addpath(genpath(pwd));
if ispc
    p = strsplit(p,';');
else
    p = strsplit(p,':');
end
job = createJob(cluster);
job.AutoAttachFiles = false;
job.AdditionalPaths = p;
end
