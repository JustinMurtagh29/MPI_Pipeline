
param = load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
param = param.p;
dataDir = fullfile(param.saveFolder, 'chiasmataSplitting/20171009T193744-kmb-on-axons-6c/outputs');
m = load(fullfile(dataDir,'20171023T115300_results.mat'));
axons = m.axons(m.indBigAxons);

batchBoundaries = round(linspace(1, numel(axons)+1, 101));

for i=1:length(batchBoundaries)-1
    batchID{i,1} = {i};
end

sharedInputs = {param,axons,batchBoundaries};
sharedInputsLocation = 2:4;
cluster = Cluster.getCluster('-p -300','-tc 20','-l h_vmem=64G','-l s_rt=3:28:00','-l h_rt=3:29:00');
job = Cluster.startJob(@connectEM.getAxonPathLength, batchID, 'name', 'pathLength', 'cluster', cluster, ...
    'sharedInputs', sharedInputs, 'sharedInputsLocation', sharedInputsLocation)

%%
%----------------------------------------------------------------

elements = dir('E:\Data_CS\Datensaetze\AxonEnding\PathLength\');
elements = elements(3:end-1);
Distance = 0;
for i=1:length(elements)
    load(strcat(elements(i).folder,'\',elements(i).name))
    Distance = Distance + sum(Distances);
end

probability =[];
for i=1:length(elements)
    load(strcat(elements(i).folder,'\',elements(i).name))
    probability = [probability; Distances/Distance];
end

IDs = [1:length(probability)]';
randomIDs = datasample(IDs,20,'Replace',false,'Weights',probability);

skelDir(randomIDs);
%%
%----------------------------------------------------------------
%{
function getAxonPathLength(batchID)

batch_folder = ls(strcat('/tmpscratch/kboerg/visX19/visX19_',int2str(batchID),'/*/result.mat'));
batch_folder = strsplit(batch_folder)';

for i=1:size(batch_folder,1)
    load(batch_folder{i,1})
    edges = minimalSpanningTree(output.nodes);
    clear distances
    for j=1:length(edges)
        distances(j,1) = pdist([output.nodes(edges(j,1),:);output.nodes(edges(j,2),:)]);
    end
    Distances(i,1) = sum(distances);
end

directory = strcat('/tmpscratch/scchr/AxonEndings/PathLength/,'int2str(batchID),'.mat');
if ~exist(directory)
    mkdir(directory)
end
save(directory, 'Distances')


function edges = minimalSpanningTree(com)
    if size(com,1) < 2
        edges = [];
    else
        % Minimal spanning tree
        adj = squareform(pdist(bsxfun(@times, com, [11.24 11.24 28])));
        tree = graphminspantree(sparse(adj), 'Method', 'Kruskal');
        [edges(:,1), edges(:,2)] = find(tree);
    end
end
%}
