rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

graphFile = fullfile(rootDir, 'graph.mat');
graph = load(graphFile, 'edges', 'borderIdx');

% NOTE(amotta): Use border idx zero instead of nan
graph.borderIdx(isnan(graph.borderIdx)) = 0;

graph = struct2table(graph);




graph(~graph.borderIdx, :) = [];

borderAreas = fullfile(rootDir, 'globalBorder.mat');
borderAreas = load(borderAreas, 'borderArea2');
borderAreas = borderAreas.borderArea2;

graph.borderArea = borderAreas(graph.borderIdx);
clear borderAreas;

% HACK(amotta): For some reason there exist borders, for which
% `physicalBorderArea2` is zero. This seems wrong.
%   In order not to be affected by this issue, let's set the area of these
% borders to NaN. This will result in a total axon-spine interface area of
% NaN, which we can remove by brute force later on.
%
% Corresponding issue on GitLab:
% https://gitlab.mpcdf.mpg.de/connectomics/auxiliaryMethods/issues/16
graph.borderArea(~graph.borderArea) = nan;
<<<<<<< HEAD
graph(:, {'prob', 'borderIdx'}) = [];
=======
graph(:, {'borderIdx'}) = [];
>>>>>>> a239cf1c32a6935ec48f6d522d4c6a5bed9838a5








save('/gaba/u/yyener/astrocyte/synapses/graph.mat', 'graph')
