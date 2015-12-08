% Load center of masses of segments and graph representation
load('/gaba/u/mberning/results/pipeline/20151111T183414/allParameter.mat');
load([p.saveFolder 'graph.mat']);
load([p.saveFolder 'CoM.mat']);

% Construct graph (now saved in graph.mat)
%edges = [graph.edges; graph.edges(:,[2 1])];
%prob = [graph.prob; graph.prob];
%graph.adj = accumarray(edges, prob, [], @max, 0, 1);
%save([p.saveFolder 'graph.mat'], 'graph', '-v7.3');

% Get seeds from skeleton for now
skel = parseNml('/gaba/u/mberning/3seeds.nml');
segIds = cell(size(skel));
for i=1:length(skel)
    pos = skel{i}.nodes(1:3) + [1 1 1];
    segIds{i} = readKnossosRoi(p.seg.root, p.seg.prefix, [pos; pos]', 'uint32', '', 'raw');
end

% Agglomerate segments and save result to skeleton
[collectedIds, probabilities] = agglomerateSG2(graph, segIds, 100);
skel = writeSkeletonEdges2(graph, com, collectedIds, probabilities);
writeNml('/gaba/u/mberning/3seedsGrown.nml', skel);

