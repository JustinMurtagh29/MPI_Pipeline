function [ tInd, nnPrec ] = calculateUpperThreshold( minPrec, synT, ...
    thresholds, scores, y_test, c_r, recLT, nnRecLT )
%CALCULATEUPPERTHRESHOLD Calculate the upper threshold based on a lower
%bound on neuron-to-neuron precision.
% INPUT minPrec: double
%           Lower bound on neuron-to-neuron precision.
%       synT: int
%           Threshold on the number of synapse per neuron pair to be
%       	considered as connected.
%       thresholds: [Nx1] double
%           Array of double containing the lower threshold in
%           the first entry and all thresholds above the lower thresholds
%           to test.
%       scores: [Nx1] double
%           Synapse detection scores on the test set.
%       y_test: [Nx1] logical
%       	Ground truth labels of test set.
%       c_r: double
%           Neuron-to-neuron connectivity ratio.
%       recLT: double
%           Synapse classifier recall for lower threshold.
%       nnRecLT: double
%           Neuron-to-neuron recall for lower threshold.
% OUTPUT tInd: int
%           Index of upper threshold in thresholds.
%        nnPrec: [Nx1] double
%           Estimated neuron-to-neuron precision all thresholds.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

lowerThres = thresholds(1);
nnPrec = zeros(length(thresholds),1);
for i = 1:length(thresholds)
    
    %calculate precision if only points above upper threshold can be fps
    TP = sum(y_test(scores >= lowerThres));
    FP = sum(~y_test(scores >= thresholds(i)));
    precision = TP/(TP + FP);
    
    %calculate nn precision
    nnFPRFA = SynEM.ErrorEstimates.nnFPR(recLT, precision,c_r, synT);
    nnPrec(i) = SynEM.ErrorEstimates.nnPrecision(nnRecLT,c_r,nnFPRFA);
end
tInd = find(nnPrec >= minPrec,1,'first');
end

