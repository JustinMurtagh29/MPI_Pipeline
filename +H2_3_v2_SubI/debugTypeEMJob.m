% restart tasks for which segmentAgglomerateFeatures.mat did not get created
% due to job failure that does not show up in job object

load('/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet/allParameterWithSynapses.mat')

cubeSize = size(p.local);
% check if all raw files exist
tic;
mainDir = fullfile(p.saveFolder, 'local');
xFiles = dir(fullfile(mainDir,'x*'));
idxError = [];
for i=1:numel(xFiles)
    thisXDir = fullfile(mainDir,xFiles(i).name);
    thisFile = fullfile(thisXDir,'segmentAgglomerateFeatures.mat');
    if ~exist(thisFile,'file')
       idxError = cat(1, idxError,i);
    end
end
toc;
Util.log(['Found ' num2str(numel(idxError)) ' tasks failed.'])

cluster = Cluster.config('memory', 48, 'time', '24:00:00', 'priority',100);
job = jobHuman; % saved in matlab session
% restart failed tasks
inputCell = {job.Tasks(idxError).InputArguments};
% Assumes all task in job execute the same function
functionH = job.Tasks(1).Function;
jobRestarted = Cluster.startJob(functionH, inputCell, 'name', 'restartedTasks', 'cluster', cluster);

