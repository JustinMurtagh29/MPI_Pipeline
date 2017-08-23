%% script to evaluate the FP synapse exclusion heuristics
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

info = Util.runInfo(false);
if ~ispc % assume on gaba
    pipelineFolder = ['/gaba/u/mberning/results/pipeline/' ...
        '20170217_ROI/'];
else % laptop
    pipelineFolder = 'E:\workspace\data\backed\20170217_ROI';
end

%% load data

Util.log('Loading SVG data.');

m = load(fullfile(pipelineFolder, 'globalEdges.mat'), 'edges');
edges = m.edges;
m = load(fullfile(pipelineFolder, 'globalSynScores.mat'));
synScores = nan(size(edges));
synScores(m.edgeIdx, :) = m.synScores;
synT = -1.67;
synEdgeIdx = max(synScores, [], 2) > -1.67;
m = load(fullfile(pipelineFolder, 'graph.mat'), 'prob', 'borderIdx');
prob = m.prob(~isnan(m.borderIdx));
m = load(fullfile(pipelineFolder, 'heuristicResult.mat'), 'myelinScore');
myelinScore = m.myelinScore;
m = load(fullfile(pipelineFolder, 'globalBorder.mat'), 'borderCoM');
borderCom = m.borderCoM;

%% get the indices of excluded predicted synapses

% prediction only based on threshold
synIdx_all = L4.Synapses.getSynapsePredictions(synScores, synT);

% prediction with exclusion of synapses with high scores for both
% interface directions
synIdx_excDir = L4.Synapses.getSynapsePredictions(synScores, synT, true);

% predicions with exclusion based on prob
synIdx_excProb = L4.Synapses.getSynapsePredictions(synScores, synT, ...
    false, prob);

% prediction with exclusion based on myelin
synIdx_excMyl = L4.Synapses.getSynapsePredictions(synScores, synT, false, ...
    [], [], edges, myelinScore);

idx_excDir = find(synIdx_all & ~synIdx_excDir);
idx_excProb = find(synIdx_all & ~synIdx_excProb);
idx_excMyl = find(synIdx_all & ~synIdx_excMyl);
n1 = length(idx_excDir);
n2 = length(idx_excProb);
n3 = length(idx_excMyl);
n = n1 + n2 + n3;

Util.log(['Discarded synapse prediction: %d based on direction, %d ' ...
    'based on prob, %d based on myelin.'], n1, n2, n3);

%% sample 60 excluded synapse prediction (20 from each criterion)

N = 33;
skel = L4.Util.getSkel();

rng('shuffle');
ridx = randi(length(idx_excDir), N, 1);
for i = 1:N
    eIdx = idx_excDir(ridx(i));
    skel = skel.addNodesAsTrees(borderCom(eIdx, :), 'DirExclusion', ...
        {sprintf('edgeIdx_%d_syn_%.3f_prob_%.3f_myl_%.3f', eIdx, ...
        max(synScores(eIdx,:)), prob(eIdx), ...
        max(myelinScore(edges(eIdx,:))))});
end

rng('shuffle');
ridx = randi(length(idx_excMyl), N, 1);
for i = 1:N
    eIdx = idx_excMyl(ridx(i));
    skel = skel.addNodesAsTrees(borderCom(eIdx, :), 'MylExclusion', ...
        {sprintf('edgeIdx_%d_syn_%.3f_prob_%.3f_myl_%.3f', eIdx, ...
        max(synScores(eIdx,:)), prob(eIdx), ...
        max(myelinScore(edges(eIdx,:))))});
end

rng('shuffle');
ridx = randi(length(idx_excProb), N, 1);
for i = 1:N
    eIdx = idx_excProb(ridx(i));
    skel = skel.addNodesAsTrees(borderCom(eIdx, :), 'ProbExclusion', ...
        {sprintf('edgeIdx_%d_syn_%.3f_prob_%.3f_myl_%.3f', eIdx, ...
        max(synScores(eIdx,:)), prob(eIdx), ...
        max(myelinScore(edges(eIdx,:))))});
end

skel.write('SynapseExclusionSamples.nml');
