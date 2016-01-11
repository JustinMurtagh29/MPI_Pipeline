% Load center of masses of segments and graph representation
load('/gaba/u/mberning/results/pipeline/20151111T183414/allParameter.mat');
load([p.saveFolder 'graph.mat']);
load([p.saveFolder 'CoM.mat']);
load([p.saveFolder 'agglomeration/nucleiVesselBorder.mat']);

% Construct graph (now saved in graph.mat)
% edges = [graph.edges; graph.edges(:,[2 1])];
% prob = [graph.prob; graph.prob];
% graph.adj = accumarray(edges, prob, [], @max, 0, 1);
%[graph.neighbours, neighbourIdx] = Graph.edges2Neighbors(graph.edges);
%for i=1:length(neighbourIdx)
%    graph.neighProb{i} = graph.prob(neighbourIdx{i});
%end
% clear neighbourIdx;

% Generate cost matrix for skeleton querrying
% edges = [graph.edges; graph.edges(:,[2 1])];
% temp = -log(graph.prob);
% cost = [temp; temp];
% graph.costM = accumarray(edges, cost, [], @min, 0, 1);
% save([p.saveFolder 'graph.mat'], 'graph', '-v7.3');

% Get seeds & GT from MH skeletons for now
folder = '/gaba/u/mberning/skeletons/source/';
files = dir([folder '*.nml']);
tr = 1;
% Read skeletons
for i=1:length(files)
    skel{i} = parseNml([folder files(i).name]);
    for j=1:length(skel{i})
        nodes{tr} = bsxfun(@plus, skel{i}{j}.nodes(:,1:3), [1 1 1]);
        tr = tr + 1;
    end
end
segIds = cell(size(nodes));
% Get segmentation ID for each node
for i=1:length(nodes)
    for j=1:size(nodes{i},1)
        pos = nodes{i}(j,:);
        segIds{i}(j) = readKnossosRoi(p.seg.root, p.seg.prefix, [pos; pos]', 'uint32', '', 'raw');
    end
end
segIds = cellfun(@(x)x(x~=0), segIds, 'UniformOutput', false);
% Get intersection of axons with border IDs
borderId = cat(1,agglo.borderMerged{:});
seeds = cellfun(@(x)intersect(x, borderId), segIds, 'UniformOutput', false);
% Fix seeds manually for now
seeds{1}(2,1) = 8787098;
seeds{2}(2,1) = 1456992;
seeds{4}(3,1) = 7320841;

% Agglomerate segments and save result to skeleton
for i=1:length(seeds)
    [collectedIds{i}, probabilities{i}, querried{i}, merger{i}, querriedEdges{i}, mergerList{i}] = agglomerateSG3(graph, com, seeds{i}, segIds{i});
    skelToWrite = writeSkeletonEdges3(graph, com, collectedIds{i}, probabilities{i}, querried{i}, merger{i}, querriedEdges{i}, skel{1}{i});
    writeNml(['/gaba/u/mberning/skeletons/axonsTest' num2str(i, '%.3i') '.nml'], skelToWrite);
end

% Visualize querried edges (analog to B4B game 1)
for i=1:length(seeds)
    for j=1:length(seeds{i})
        queryIds = querriedEdges{i}{j}(:,1);
        for k=1:length(queryIds)
            queryIdx = find(collectedIds{i}{j} == queryIds(k));
            aggloIds = collectedIds{i}{j}(1:queryIdx);
            % TO DO: Add handling of merger here & in function
            mergerIds = [];
            outputFile = ['/gaba/u/mberning/skeletons/movies/queries' num2str(i, '%.2i') '_' num2str(j, '%.2i') '_' num2str(k, '%.2i')];
            arbitraryResliceAgglo(p, graph, com, aggloIds, queryIds(k), mergerIds, outputFile);
        end
    end
end


