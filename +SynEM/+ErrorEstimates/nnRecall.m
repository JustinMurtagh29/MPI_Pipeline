function p = nnRecall( synR, synT, empDist )
%NTONRECALL Calculate the neuron-to-neuron recall from synapse recall.
% The n-n recall is given by the probability of finding at least synT
% synapses between these neurons given that each synapse is detected
% independently by a retrieval probability. An empirical distribution over
% the number of synapses per neuron pair is used.
% INPUT synR: [Nx1] double
%           Probabilities of detecting a synapse (e.g. the recall of the
%           synapse classifier on the test set).
%       synT: int
%       	Number of synapses to accept a neuron-to-neuron connection.
%       empDist: (Optional) [1xN] double
%           Vector specifying the distribution of number of synapse between
%           a connected neuron pair. The i-th entry specifies the
%           probability that a connected pair of neurons has i synapses.
%           (Default: Empirical distribution from Feldmeyer paper).
% OUTPUT p: [Nx1] array of double containing the neuron-to-neuron recall
%           for each retrieval probability.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('empDist','var') || isempty(empDist)
    empDist = [0 2 5 2 2]./11;
else %normalize to make sure
    empDist = empDist./sum(empDist);
end

p = zeros(length(synR),1);
for kappa = 1:length(empDist)
    p = p + binocdf(synT - 1, kappa, synR,'upper').*empDist(kappa);
end


end

