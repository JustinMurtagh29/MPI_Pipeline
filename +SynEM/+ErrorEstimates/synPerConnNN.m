function [ empDist, n_syn ] = synPerConnNN( type, normalize )
%SYNPERCONNNN Returns empirical distributions for the number of synapses
%between connected neurons for different neuron types/layers of mouse/rat
%somatosensory cortex.
% INPUT type: string
%           Type of the distribution. Options are
%               'combined'
%               'L4exc-L4exc'
%       normalize: (Optional) logical
%           Flag indicating whether that the relative frequencies of
%           synapses are used.
% OUTPUT empDist: [Nx1] double
%           Cumulative or relative frequency distribution where the i-th
%           entry corresponds to i synapses.
%        n_syn: double
%           Average number of synapses.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('normalize', 'var') || isempty(normalize)
    normalize = true;
end

switch type
    case 'combined'
        empDist = [1 4 13 11 19 5 3 1]; %combined for all layers
    case 'L4exc-L4exc'
        empDist = [0 2 5 2 2]; %this are L4-L4 connections
    case 'L4inh-L4'
        %Koelbl et al., 2015, CerebCortex
        empDist = [0 4 5 5 2 1];
    otherwise
        error('Unknown type %s.', type);
end

if normalize
    empDist = empDist./sum(empDist);
end

n_syn = sum(empDist./sum(empDist(:)).*(1:length(empDist)));

end

