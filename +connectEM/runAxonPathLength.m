%% Axons
param = load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
param = param.p;
dataDir = fullfile(param.saveFolder, 'chiasmataSplitting/20171009T193744-kmb-on-axons-6c/outputs');
m = load(fullfile(dataDir,'20171023T115300_results.mat'));
axons = m.axons(m.indBigAxons);

batchBoundaries = round(linspace(1, numel(axons)+1, 101));

for i=1:length(batchBoundaries)-1
    batchID{i,1} = {i};
end

dataDir = fullfile(param.saveFolder, 'aggloState/axonPathLength/');
if ~exist(dataDir)
    mkdir(dataDir)
end

sharedInputs = {param,axons,batchBoundaries,dataDir};
sharedInputsLocation = 2:4;
cluster = Cluster.getCluster('-p -300','-tc 20','-l h_vmem=64G','-l s_rt=3:28:00','-l h_rt=3:29:00');
job = Cluster.startJob(@connectEM.getAxonPathLength, batchID, 'name', 'pathLength', 'cluster', cluster, ...
    'sharedInputs', sharedInputs, 'sharedInputsLocation', sharedInputsLocation)

%%
%----------------------------------------------------------------
aggloLengths=[];
totalPathLength = 0;
for i=1:100
    load(fullfile(param.saveFolder, strcat('aggloState/axonPathLength/batch',num2str(i),'.mat')))
    totalPathLength = totalPathLength + sum(aggloLength);
    aggloLengths = cat(1,aggloLengths,aggloLength);
end

probabilities = aggloLengths./totalPathLength;

IDs = [1:length(probabilities)]';
randomIDs = datasample(IDs,50,'Replace',false,'Weights',probabilities);

folder = fullfile(param.saveFolder, 'aggloState/axonPathLength/randomExamples');
for i=1:length(randomIDs)
    randomAgglo = axons(randomIDs(i));

    connectEM.generateSkeletonFromAggloNew(randomAgglo, {strcat('Axon',num2str(randomIDs(i)),'PathLength_',num2str(aggloLengths(randomIDs(i))))} , folder, [],[],strcat('randomAxon',num2str(randomIDs(i)),'.nml'));
end
%% Dendrites:
param = load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
param = param.p;
dataDir = fullfile(param.saveFolder, 'aggloState');
m = load(fullfile(dataDir,'dendrites_07.mat'));
dendrites = m.dendrites(m.indBigDends);

batchBoundaries = round(linspace(1, numel(dendrites)+1, 101));

for i=1:length(batchBoundaries)-1
    batchID{i,1} = {i};
end

dataDir = fullfile(param.saveFolder, 'aggloState/dendritePathLength/');
if ~exist(dataDir)
    mkdir(dataDir)
end

sharedInputs = {param,dendrites,batchBoundaries,dataDir};
sharedInputsLocation = 2:4;
cluster = Cluster.getCluster('-p -300','-tc 20','-l h_vmem=64G','-l s_rt=3:28:00','-l h_rt=3:29:00');
job = Cluster.startJob(@connectEM.getAxonPathLength, batchID, 'name', 'pathLength', 'cluster', cluster, ...
    'sharedInputs', sharedInputs, 'sharedInputsLocation', sharedInputsLocation)

for i=1:length(batchBoundaries)-1
    connectEM.getAxonPathLength(i,param,dendrites,batchBoundaries,dataDir);
end

%%
%----------------------------------------------------------------
aggloLengths=[];
totalPathLength = 0;
for i=1:100
    load(fullfile(param.saveFolder, strcat('aggloState/dendritePathLength/batch',num2str(i),'.mat')))
    totalPathLength = totalPathLength + sum(aggloLength);
    aggloLengths = cat(1,aggloLengths,aggloLength);
end

probabilities = aggloLengths./totalPathLength;

IDs = [1:length(probabilities)]';
randomIDs = datasample(IDs,50,'Replace',false,'Weights',probabilities);

folder = fullfile(param.saveFolder, 'aggloState/axonPathLength/randomExamples');
for i=1:length(randomIDs)
    randomAgglo = axons(randomIDs(i));

    connectEM.generateSkeletonFromAggloNew(randomAgglo, {strcat('Axon',num2str(randomIDs(i)),'PathLength_',num2str(aggloLengths(randomIDs(i))))} , folder, [],[],strcat('randomAxon',num2str(randomIDs(i)),'.nml'));
end
