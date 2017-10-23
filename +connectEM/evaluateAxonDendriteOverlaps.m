% This function takes a bunch of axon super-agglomerates and checks if the
% flight paths contained therein have significant overlaps with some of the
% dendrites agglomerates.
%
% The hope is that this helps us to identify axon fragments which ended up
% in the wrong bucket. It could be that these fragments are a major source
% of "dendrite" endings. Identifying these fragments could potentially help
% us to get rid of unnecessary queries.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>