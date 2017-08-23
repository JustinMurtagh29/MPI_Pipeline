% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
% script comparing a synapse in the original segmentation and in the
% correspondence agglomeration

%% load svg data
load('/gaba/u/bstaffle/code/workspace/allParameter20170217.mat')
m = load([p.saveFolder 'globalEdges.mat']);
edges = m.edges;
m = load([p.saveFolder 'correspondencesNew.mat']);
corrEdges = m.corrEdges;
m = load([p.saveFolder 'globalSynScores.mat']);
synScores = nan(size(edges));
synScores(m.edgeIdx, :) = m.synScores;

%% define edges
pos = [640, 1177, 385]; % to find it in wk
idsPre = [1386895, 1395212, 176709, 184835];
idsPost = [176653, 184723, 1395170, 1386874];

%% define edges

pos = [1152, 4890, 914];
idsPre = [4264728, 4271617, 4264730];
idsPost = [4271008, 4264512];

%% define edges

pos = [2456, 4221, 871];
idsPre = [2932408, 2932407, 3012096, 3012173, 2932111, 3011580];
idsPost = 2932347;

%% define edges

pos = [5225, 785, 3201];
idsPre = [14089466, 14089493, 14089242, 14095480, 14095481, 12908137, ...
    12908123, 12908123];
idsPost = [14089188, 12903552];

%% define edges

pos = [5248, 735, 3407];
idsPre = [14094729, 14099497];
idsPost = [14094810, 14099438];

%% get edges between pre and post

[~,edgeIdx] = L4.Agglo.findEdgesBetweenAgglos({idsPre, idsPost}, edges);
edgeIdx = edgeIdx{1};
thisSynScores = synScores(edgeIdx,:);

%% get agglo border mapping

db = Codat.StoreEM(p);
aggloName = 'correspondences';
aggloFolder = db.getAggloFolder(aggloName);
m = load([aggloFolder 'edges.mat'], 'borderMapping');
borderMapping = m.borderMapping;
m = load([aggloFolder 'synapseScores.mat']);
synScoresAgglo = nan(max(borderMapping), 2);
synScoresAgglo(m.edgeIdx, :) = m.scores;

%% edges in agglo

aggloEdgeIdx = unique(borderMapping(edgeIdx));
thisSynScoresAgglo = synScoresAgglo(aggloEdgeIdx,:);

%% print output

Util.log(['Score of synapses changed from %.4f in original ' ...
    'segmentation to %.4f in correspondence agglomeration'], ...
    max(thisSynScores(:)), max(thisSynScoresAgglo(:)));
