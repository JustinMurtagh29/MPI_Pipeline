which_col = 1; %the super merger
temp{which_col} = load(['/tmpscratch/kboerg/visX11/visX11_' num2str(floor(which_col/100)) '/visX11_' num2str(which_col) '/result.mat']);


outputFolder = '/tmpscratch/kboerg/';
scratchFolder = outputFolder;
skeletonFolders = {'MBKMB_L4_chiasma_axon_queries_2017_a_nmls' 'MBKMB_L4_chiasma_axon_queries_2017_b_nmls'};
skeletonFolders = cellfun(@(x)[scratchFolder x filesep], skeletonFolders, 'uni', 0);

temp2 = load('/gaba/scratch/mberning/axonQueryGeneration/beforeQueryGeneration.mat', 'axonsNew');
axons = temp2.axonsNew;
clear temp2;
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');

[ff.segIds, ff.neighbours, ff.filenames, ff.nodes, ff.startNode, ff.comments] = connectEM.lookupNmlMulti(p, skeletonFolders, false);
ff = structfun(@(x)x(cellfun(@isempty,ff.comments)), ff, 'uni', 0);

idxtemp = cellfun(@(x){find(ismember(cell2mat(ff.startNode'), [-1,-1,-1;cell2mat(x.output.position)],'rows'))},temp);
fflookup = repelem(1:length(ff.startNode),cellfun(@(x)size(x,1),ff.startNode));
p.sphereRadiusOuter = Inf; % in nm
p.sphereRadiusInner = 1000; % in nm
p.voxelSize = [11.24 11.24 28];

nodesScaled = bsxfun(@times,temp{which_col}.output.nodes,p.voxelSize);
skelsegids = connectEM.lookupSkelGT(p, struct('nodes',{{temp{which_col}.output.nodes}}));
save('/tmpscratch/kboerg/chiasmarun/initial_segids.mat','skelsegids');
clear skelsegids
save('/tmpscratch/kboerg/chiasmarun/initial.mat');

functionH = @connectEM.splitChiasmataSub;
inputCell = cellfun(@(x){x}, num2cell(1:2000), 'uni', 0);
cluster = Cluster.getCluster( ...
    '-pe openmp 1', ...
    '-p 0', ...
    '-l h_vmem=24G', ...
    '-l s_rt=23:50:00', ...
    '-l h_rt=24:00:00');
job = Cluster.startJob( functionH, inputCell, ...
    'name', 'chiasmata', ...
    'cluster', cluster);    
%% wait for job to be done
functionH = @connectEM.splitChiasmaSub2;
inputCell = cellfun(@(x){x}, num2cell(1:500), 'uni', 0);

job = Cluster.startJob( functionH, inputCell, ...
    'name', 'chiasmata', ...
    'cluster', cluster);    

thisEdgesNewMeta = [];
thisEdgesColMeta = temp{which_col}.output.edges;
listmats = dir('/tmpscratch/kboerg/chiasmarun2/resultchiasma2_*');
for idx = 1:length(listmats)
    idx
    load(['/tmpscratch/kboerg/chiasmarun2/' listmats(idx).name]);
    thisEdgesNewMeta=[thisEdgesNewMeta;thisEdgesNew];
    thisEdgesColMeta=intersect(thisEdgesColMeta, thisEdgesCol, 'rows');
end
thisEdgesFinal=[thisEdgesColMeta;thisEdgesNewMeta];
C = Graph.findConnectedComponents(thisEdgesFinal);
