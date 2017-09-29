function job = splitChiasmata()
result3 = load('/tmpscratch/kboerg/axonQueryAnalysisResult3.mat');
lookup_tasks = load('/tmpscratch/kboerg/chiasmarunAugust/lookup_tasks');
lookup_tasks = lookup_tasks.lookup_tasks;

outputFolder = '/tmpscratch/kboerg/chiasmarunAugust/';
scratchFolder = outputFolder;

% NOTE(amotta): path to directories with NML of chiasmata queries
skeletonFolders = {'nmls4'};
skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);

p = load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
p = p.p;

axons = load('/tmpscratch/kboerg/axonsBorderNew');
axons = axons.axons;

%% Start doing actual work
% TODO(amotta): Shouldn't we ignore NMLs with multiple trees?
[ff.segIds, ff.neighbours, ff.filenames, ff.nodes, ff.startNode, ff.comments, ff.time] = connectEM.lookupNmlMulti(p, skeletonFolders, false);

% NOTE(amotta): Ignore flight paths comments
ff = structfun(@(x)x(cellfun(@isempty,ff.comments)), ff, 'uni', 0);

% NOTE(amotta): I presume that the below `ids2.txt` file contains a list
% with the IDs of all created (and to be processed) tasks. First, this CSV
% file is parsed and the task IDs are extracted.
allbutfirst = @(x)x(2:end);
first = @(x)x(1);
fid = fopen('/tmpscratch/kboerg/chiasmarunAugust/ids2.txt');
M = cellfun(@(x)first(strsplit(x,',')),allbutfirst(strsplit(char(fread(fid))','\n')));
fclose(fid);

% NOTE(amotta): Try to find the task IDs in the list of parsed NML files.
% `idx_pre_pre` is a logical vector indicating which of the handed out
% tasks were found again.
%
% `lookup_tasks` is a matrix which, presumably`, contains
% * `lookup_tasks(:, 1)`, the agglomerate ID
% * `lookup_tasks(:, 2)`, the chiasma ID within the agglomerate
%
% `idx_pre` then contains the list of all agglomerates that are affected by
% the chiasmata splitting procedure.
idx_pre_pre = cellfun(@(x)max(cellfun(@length,strfind(ff.filenames,x))),M)>0;
idx_pre = unique(lookup_tasks(idx_pre_pre,1));

% NOTE(amotta): Save all variable for use in sub function
save('/tmpscratch/kboerg/chiasmarunAugust/preSplitRun');

% NOTE(amotta): One task per 50 chiasmata!
functionH = @connectEM.splitChiasmataSub;
inputCell = cellfun(@(x){x}, num2cell(1 : 50), 'uni', 0);

cluster = Cluster.getCluster( ...
    '-pe openmp 1', ...
    '-p 0', ...
    '-l h_vmem=24G', ...
    '-l s_rt=23:50:00', ...
    '-l h_rt=24:00:00');
job = Cluster.startJob( functionH, inputCell, ...
    'name', 'chiasmata2', ...
    'cluster', cluster);
end
function dummy()
end