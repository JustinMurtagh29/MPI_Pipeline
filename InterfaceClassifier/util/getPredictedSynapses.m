function [centroidList] = getPredictedSynapses(s,cubes)
%GETPREDICTEDSYNAPSES Get the centroids of all predicted synapses with score bigger than s.
% INPUT s: Score threshold. All interfaces with score bigger than s are considered synaptic.
%       cubes: p.local(:) from allParameters
% OUTPUT centroidList: Cell-array with the centroids of all synapses.
centroidList = cell(0); 
for i = 1:length(cubes)
    m = matfile(cubes(i).synapseFile);
    scores = m.scores;
    centroids_cube = m.centroidList;
    synapses = scores > s;
    synapses_single = synapses(1:end/2) | synapses(end/2 + 1:end);
    bbox = cubes(i).bboxSmall(:,1)';
    synapseCentroids = centroids_cube(synapses_single);
    synapseCentroids = cellfun(@(x) x + bbox,synapseCentroids,'UniformOutput',false);
    centroidList = cat(1,centroidList,synapseCentroids);
end 
end
