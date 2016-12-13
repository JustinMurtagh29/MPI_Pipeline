function [tInd, recP] = calculateLowerThreshold( minRec, synT, ...
    recallValues, empDist )
%CALCULATELOWERTHRESHOLD Calculate the lower score threshold based on a
%lower bound on neuron-to-neuron recall.
% INPUT minRec: double
%           The lower bound on the neuron-to-neuron recall.
%       synT: Double
%           Threshold on the number of synapses per neuron pair which
%       	need to be detected (i.e. how many synapses between
%       	to connections we need in order to accept these neurons as
%           connected).
%       recallValues: [Nx1] double
%           Recall values of the classifier.
%       empDist: (Optional) [1xN] double
%           Vector specifying the distribution of number of synapse between
%           a connected neuron pair. The i-th entry specifies the
%           probability that a connected pair of neurons has i synapses.
%           (Default: Empirical distribution from Feldmeyer paper).
% OUTPUT tInd: The lower threshold index as the index of the
%        	corresponding recall value.
%        recP: The resulting neuron-to-neuron recall for each synapse
%        	recall.
%
% see also nnRecall.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('empDist','var')
    empDist = [];
end

pValues = zeros(length(recallValues),1);
for i = 1:length(recallValues)
    pValues(i) = SynEM.ErrorEstimates.nnRecall(recallValues(i), ...
        synT, empDist);
end
recP = pValues;
pValues(pValues < minRec) = inf;
[~,tInd] = min(pValues);

end