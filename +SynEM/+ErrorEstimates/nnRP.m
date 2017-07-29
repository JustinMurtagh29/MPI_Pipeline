function [ nnRP ] = nnRP( rp, c_r, n_syn, synT, empDist )
%NNRP Calculate neuron-to-neuron RP from synapse RP.
% INPUT rp: [Nx2] double
%           Recall-precision value pairs for synapse detection.
%       c_r: double
%           Neuron-to-neuron connectivity ratio.
%       n_syn: double
%           Average number of synapses between two connected neurons.
%       synT: (Optional) int
%           Number of synapses that is required in order to accept a
%           neuron-to-neuron connection as detected.
%           (Default: 1)
%       empDist: (Optional) [1xN] double
%           Vector specifying the distribution of number of synapse between
%           a connected neuron pair. The i-th entry specifies the
%           probability that a connected pair of neurons has i synapses.
%           (Default: Empirical distribution from Feldmeyer paper).
% OUTPUT nnRP: [Nx2] double
%           Neuron-to-neuron recall-precision value pairs for each synapse
%           rp value pair.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

synR = rp(:,1);
synP = rp(:,2);
nnR = SynEM.ErrorEstimates.nnRecall(synR, synT, empDist);
nnP = SynEM.ErrorEstimates.nnPrecision(synR, synP, c_r, synT, ...
    n_syn, nnR);
nnRP = [nnR, nnP];

end

