function synWKMapping( edges, synScores, t )
%SYNWKMAPPING Synapse visualization mapping for webknossos.
% INPUT edges: [Nx2] array of integer specifying the global edge list.
%       synScores: [Nx2] array of double specifying the syn scores for each
%           edge.
%       t: Threshold for synapse detection.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

synScores = max(synScores, [], 2);
syn = edges(synScores > t, :);
other = edges(synScores <= t, :);
other = setdiff(unique(other(:)), syn(:));
components = num2cell(syn,2);
components = cat(1,components, {[0;other]});
components = cellfun(@double,components,'UniformOutput',false);
wk.makeWKMapping(components, 'synapses');

end