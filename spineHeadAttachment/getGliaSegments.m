function gliaSegments = getGliaSegments(graph)
load('/gaba/u/mberning/results/pipeline/20151111T183414/gliaSegments.mat');
gliaSegments = getXSegments(graph, glia, 0.96, Inf, NaN, Inf, []);