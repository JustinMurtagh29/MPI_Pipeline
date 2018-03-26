% This script exports axon-dendrite pairs that are coupled by exactly N
% synapses according to the connectome. These candidates are then inspected
% in webKNOSSOS and manually labelled as true or false.
%
% In the "true" case, the pre- and postsynaptic volumes are manually
% reconstructed by picking up the corresponding segments. This will give a
% good estimate of the true contact area.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;
