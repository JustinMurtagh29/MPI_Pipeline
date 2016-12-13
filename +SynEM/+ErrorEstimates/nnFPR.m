function p_fpr = nnFPR( synR, synP, c_r, synT, n_syn )
%NNFPR Probability that not-connected neuron pair gets connected.
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
% OUTPUT p_fpr: [Nx1] double
%           Probability of making a false positive neuron-to-neuron
%           connection for the corresponding row in synR and synP.
%
% NOTE synR and synP must have the same length.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('synT','var') || isempty(synT)
    synT = 1;
end
if ~exist('n_syn','var') || isempty(n_syn)
    empDist = [0 2 5 2 2]./11;
    n_syn = sum(empDist.*(1:length(empDist)));
end
if synP == 0
    p_fpr = 1;
    return
end
lambda = (1-synP).*synR./synP.*c_r.*n_syn;
p_fpr = poisscdf(synT - 1,lambda,'upper');

end

