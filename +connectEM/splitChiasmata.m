function job = splitChiasmata()
result3 = load('/tmpscratch/kboerg/axonQueryAnalysisResult3.mat');
lookup_tasks = load('/tmpscratch/kboerg/chiasmarunAugust/lookup_tasks');
lookup_tasks = lookup_tasks.lookup_tasks;

outputFolder = '/tmpscratch/kboerg/chiasmarunAugust/';
scratchFolder = outputFolder;
skeletonFolders = {'nmls4'};
skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);


p = load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
p = p.p;
axons = load('/tmpscratch/kboerg/axonsBorderNew');
axons = axons.axons;
[ff.segIds, ff.neighbours, ff.filenames, ff.nodes, ff.startNode, ff.comments, ff.time] = connectEM.lookupNmlMulti(p, skeletonFolders, false);
ff = structfun(@(x)x(cellfun(@isempty,ff.comments)), ff, 'uni', 0);
allbutfirst = @(x)x(2:end);
first = @(x)x(1);
 fid = fopen('/tmpscratch/kboerg/chiasmarunAugust/ids2.txt');
 M = cellfun(@(x)first(strsplit(x,',')),allbutfirst(strsplit(char(fread(fid))','\n')));
 fclose(fid);
 idx_pre_pre = cellfun(@(x)max(cellfun(@length,strfind(ff.filenames,x))),M)>0;
 idx_pre = unique(lookup_tasks(idx_pre_pre,1));

save('/tmpscratch/kboerg/chiasmarunAugust/preSplitRun');
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