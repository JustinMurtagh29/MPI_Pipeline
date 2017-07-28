function [ empDist, n_syn, c_r ] = synPerConnNN( type, normalize )
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
%        c_r: double
%           Connectivity ratio for specified cell pairs (if available)
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('normalize', 'var') || isempty(normalize)
    normalize = true;
end

switch type
    case 'exc_combined'
        empDist = [1 4 13 11 19 5 3 1]; %combined for all layers
        c_r = 0.2;
    case 'inh_combined'
        empDist = [0 0 0 0 0 1 0 0 0 0];
        c_r = 0.6;
    case 'L4exc-L4exc'
        %Feldmeyer et al., 1999
        empDist = [0 2 5 2 2]; %this are L4-L4 connections
        c_r = 0.2;
    case {'L4inh-L4', 'Koelbl2015'}
        %Koelbl et al., 2015, CerebCortex
%         empDist = [0 4 5 5 2 1]; %all (see paper)
        empDist = [0 2 2 4 1 1]; %barrel confined
        c_r = 0.67;
    case {'L23inh-L23pyr', 'Hoffmann2015'}
        %empDist ranging from 3-10 but distr not in paper
        %use something rather flat that matches the statistics
        empDist = [0 0 1 3 4 2 2 2 2 1]./17;
        c_r = 0.18; %ranging from 0.16-0.21 (see paper)
    otherwise
        error('Unknown type %s.', type);
end

if normalize
    empDist = empDist./sum(empDist);
end

if ~exist('c_r', 'var')
    c_r = [];
end

n_syn = sum(empDist./sum(empDist(:)).*(1:length(empDist)));

end

