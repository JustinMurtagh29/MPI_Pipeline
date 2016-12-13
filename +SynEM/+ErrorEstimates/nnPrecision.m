function [ nnP ] = nnPrecision( synR, synP, c_r, synT, n_syn, nnR )
%NNPRECISION Neuron-to-neuron precision.
% INPUT synR: [Nx1] double
%           Synapse recall value.
%       synP: [Nx1] double
%           Synapse precision values.
%       c_r: double
%           Neuron-to-neuron connectivity ratio.
%       synT: (Optional) int
%           Number of synapses that is required in order to accept a
%           neuron-to-neuron connection as detected.
%           (Default: 1)
%       n_syn: (Optional) double
%           Average synapse number per connection.
%           (Default: Average from Feldmeyer distribution).
%       nnR: [Nx1] double
%           Neuron-to-neuron recall values.
%           (see output of SynEM.ErrorEstimates.nnRecall)
% OUTPUT nnP: [Nx1] double of neuron-to-neuron precision rates for the
%           respective rows in synR, synP and nnR.
%
% NOTE synR, synP and nnR must have the same length.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

nnFPR = SynEM.ErrorEstimates.nnFPR( synR, synP, c_r, synT, n_syn );
nnP = SynEM.ErrorEstimates.nnPrecisionFromFPR(nnR, c_r, nnFPR);

end

