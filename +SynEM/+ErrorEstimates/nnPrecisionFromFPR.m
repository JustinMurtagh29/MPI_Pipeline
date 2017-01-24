function [ nnPrec ] = nnPrecisionFromFPR( nnRec, c_r, nnFPR )
%NNPRECISIONFROMFPR Neuron-to-neuron precision from nnFPR.
% INPUT nnRec: [Nx1] double
%           Neuron-to-neuron recall values.
%           (see output of SynEM.ErrorEstimates.nnRecall)
%       c_r: double
%           Neuron-to-neuron connectivity ratio.
%       nnFPR: [Nx1] double
%           Neuron-to-neuron false positive probabilities.
%       	(see output of SynEM.ErrorEstimates.nnFPR).
% OUTPUT nnPrec: Neuron-to-neuron precision.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

nnPrec = nnRec./(nnRec + (1 - c_r)/c_r * nnFPR);

end