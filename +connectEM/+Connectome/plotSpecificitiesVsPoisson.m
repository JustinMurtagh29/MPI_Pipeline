% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

connFile = fullfile( ...
    rootDir, 'connectomeState', ...
    'connectome_axons_18_a_with_den_meta.mat');

minSynPre = 10;
info = Util.runInfo();

%% loading data
conn = load(connFile);

%% calculate target-class probability over synapses
[targetClasses, ~, targetClassSyns] = unique(conn.denMeta.targetClass);
targetClassSyns = accumarray(targetClassSyns, conn.denMeta.synCount);
targetClassProbs = targetClassSyns ./ sum(targetClassSyns);

%% restrict to axons with enough synapses
axonIds = find(conn.axonMeta.synCount >= minSynPre);
[synCounts, ~, synCountAxons] = unique(conn.axonMeta.synCount(axonIds));
synCountAxons = accumarray(synCountAxons, 1);

className = 'ApicalDendrite';
classProb = targetClassProbs(targetClasses == className);

poiss = table;
poiss.prob = cell2mat(arrayfun(@(nSyn, nAxons) ...
    nAxons * poisspdf((0:nSyn)', nSyn * classProb), ...
    synCounts, synCountAxons, 'UniformOutput', false));
poiss.spec = cell2mat(arrayfun( ...
    @(nSyn) (0:nSyn)' ./ nSyn, ...
    synCounts, 'UniformOutput', false));


binEdges = linspace(0, 1, 21);
poiss.binId = discretize(poiss.spec, binEdges);
binCount = accumarray(poiss.binId, poiss.prob);

figure;
histogram('BinEdges', binEdges, 'BinCounts', binCount)