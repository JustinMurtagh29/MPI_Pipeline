function [ centroids ] = calculateSynapseCentroids( data, interfaceLabels )
%CALCULATESYNAPSECENTROIDS Calculate centroids of synapses in data.
% INPUT data: Interface classifier data struct.
%       interfaceLabels: Label vector for interfaces in data.
% OUTPUT centroids: n x 3 array with centroids coordinates
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

synapses = find(interfaceLabels < 3);
centroids = zeros(length(synapses),3);
for i = 1:length(synapses)
    [x,y,z] = ind2sub(size(data.raw),data.interfaceSurfaceList{synapses(i)});
    centroids(i,:) = [mean(x),mean(y),mean(z)];
end