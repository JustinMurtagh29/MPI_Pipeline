% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>

clear;
info = Util.runInfo();

rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
synapseFile = fullfile(rootDir, 'connectomeState/SynapseAgglos_v2.mat');
connectomeFile = fullfile(rootDir, 'connectomeState/connectome.mat');
outFile = '/home/amotta/Desktop/connectomeAugmented.mat';

%% load input data
synData = load(synapseFile);
data = load(connectomeFile);

%% add number of synapses per edge
data.connectome.synCount = cellfun( ...
    @numel, data.connectome.synIdx);

%% per axon statistics
data.axonMeta = table;
data.axonMeta.id = reshape(1:numel(data.axons), [], 1);
data.axonMeta.synCount = accumarray( ...
    data.connectome.edges(:, 1), ...
    cellfun(@numel, data.connectome.synIdx), ...
   [numel(data.axons), 1], @sum, 0);
data.axonMeta.spineSynCount = accumarray( ...
    data.connectome.edges(:, 1), cellfun(@(ids) ...
    sum(synData.isSpineSyn(ids)), data.connectome.synIdx), ...
   [numel(data.axons), 1], @sum, 0);

%% target classes
data.targetClass = repelem( ...
    {'OtherDendrite'}, numel(data.dendrites), 1);
data.targetClass(data.idxAD) = {'ApicalDendrite'};
data.targetClass(data.idxSD) = {'SmoothDendrite'};
data.targetClass(data.idxSoma) = {'Soma'};
data.targetClass = categorical(data.targetClass);

%% save output
data.info = info;
Util.saveStruct(outFile, data);