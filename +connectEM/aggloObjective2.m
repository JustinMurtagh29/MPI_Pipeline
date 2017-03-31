function [value, goodWeight, badWeight, nrEdgesGood, nrEdgesBad] = aggloObjective2( graph, cc )
% Calculate probability within and between connected components

% Design parameter (higher means more oversegmentation)
beta = 0.3;

% Get all edges & corresponding prob & weights in component
idx = all(ismember(graph.edges, cc),2);
theseProbabilities =  max(graph.prob(idx), 0.01);
theseWeights = log(theseProbabilities) + log((1-beta)./beta);

% Sum over positive and negative edges respectively
goodWeight = sum(theseWeights(theseWeights > 0));
badWeight = sum(theseWeights(theseWeights < 0));
value = goodWeight + badWeight;

% Number of edges with good or bad weight
nrEdgesGood = sum(theseWeights > 0);
nrEdgesBad = sum(theseWeights < 0);


end
