%% script to generate the connectome
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%% load synapse predictions

m = load('E:\workspace\data\backed\20170217_ROI\globalEdges.mat', 'edges');
edges = m.edges;

m = load('E:\workspace\data\backed\20170217_ROI\globalSynScores.mat');
synScores = nan(size(edges));
synScores(m.edgeIdx, :) = m.synScores;

%% get agglos

m = load('E:\workspace\data\backed\L4\FocusEM\postQueryAgglomerates.mat');
% m.isAbove5um(m.toIgnoreIdx) = false;
% axonsPostQuery = m.axonsPostQuery(m.isAbove5um);
axonsPostQuery = m.axonsPostQuery;
m = load('E:\workspace\data\backed\L4\FocusEM\dendritesPostQuery.mat');
% dendritesPostQuery = m.dendritesPostQuery(m.isAbove5um);
dendritesPostQuery = m.dendritesPostQuery;

%% get edges between agglos

agglos = cat(1, axonsPostQuery, dendritesPostQuery);
[c, edgeIdx] = L4.Agglo.findEdgesBetweenAgglos(agglos, edges);

% number of synapses for corresponding edges
numSyn = cellfun(@(x)sum(any(synScores(x, :) > -1.67, 2)), edgeIdx);

%% write to connectome

% detele edges between axons and dendrites
num_ax = length(axonsPostQuery);
num_de = length(dendritesPostQuery);
toDel = all(c <= num_ax, 2) | all (c > num_ax, 2);
c(toDel,:) = []; %first row ax, second dendr
numSyn(toDel) = [];
C = zeros(num_ax, num_de);
idx = sub2ind(size(C), c(:,1), c(:,2) - num_ax);
C(idx) = numSyn;
