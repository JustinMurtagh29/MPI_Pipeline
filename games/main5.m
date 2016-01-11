% Load center of masses of segments and graph representation
load('/gaba/u/mberning/results/pipeline/20151111T183414/allParameter.mat');
graph = load([p.saveFolder 'graph.mat']);
boutons = load([p.saveFolder 'Boutons.mat']);
load([p.saveFolder 'CoM.mat']);
load([p.saveFolder 'agglomeration/nucleiVesselBorder.mat']);

% Get seeds & GT from MH skeletons for now
folder = '/gaba/u/mberning/skeletons/source/biggerSet/';
files = dir([folder '*.nml']);
for i=1:length(files)
    skel{i} = parseNml([folder files(i).name]);
    if length(skel{i}) > 1
        error('More than one tree in nml file');
    end
    if unique(skel{i}{1}.edges(:))' ~= 1:size(skel{i}{1}.nodes,1)
        error('Node missing in edge list');
    end
    nodeDegree{i} = histc(skel{i}{1}.edges(:),unique(skel{i}{1}.edges(:)));
    nodes{i} = bsxfun(@plus, skel{i}{1}.nodes(:,1:3), [1 1 1]);
end

% QUICK FIX: use only skeletons with 5 or less nodes of degree one
excludeIdx = cellfun(@(x)sum(x==1), nodeDegree) > 5;
skel(excludeIdx) = [];
nodes(excludeIdx) = [];
nodeDegree(excludeIdx) = [];
files(excludeIdx) = [];
% QUICK FIX: remove skeletons with less than 20 nodes
excludeIdx = cellfun(@(x)size(x,1), nodeDegree) < 20;
skel(excludeIdx) = [];
nodes(excludeIdx) = [];
nodeDegree(excludeIdx) = [];
files(excludeIdx) = [];

segIds = cell(size(nodes));
for i=1:length(nodes)
    for j=1:size(nodes{i},1)
        pos = nodes{i}(j,:);
        segIds{i}(j) = readKnossosRoi(p.seg.root, p.seg.prefix, [pos; pos]', 'uint32', '', 'raw');
    end
end

% Restrict segmentation IDs to those not equal to 0
segIdsR = cellfun(@(x)x(x~=0), segIds, 'UniformOutput', false);
% Get intersection of axons with border IDs
borderId = cat(1,agglo.borderMerged{:});
seeds = cellfun(@(x)intersect(x, borderId), segIdsR, 'UniformOutput', false);

% QUICK FIX: remove skeletons where no seed was found
excludeIdx = cellfun(@(x)size(x,1), seeds) == 0;
skel(excludeIdx) = [];
nodes(excludeIdx) = [];
nodeDegree(excludeIdx) = [];
files(excludeIdx) = [];
seeds(excludeIdx) = [];

% Chose one random seed for each axon
nrSeeds = cellfun(@length, seeds);
for i=1:length(seeds)
    randIdx = randi(nrSeeds(i), 1);
    seed{i} = seeds{i}(randIdx);
end

% Agglomerate segments and save result to skeleton
matlabpool 8;
t_prob = .75:0.01:0.99;
t_neigh = 20:5:40;
tic;
for t1=1:length(t_prob)
    for t2=1:length(t_neigh)
        parfor i=1:length(seed)
            [collectedIds{t1,t2}{i}, probabilities{t1,t2}{i}, mergerList{t1,t2}{i}, queryId{t1,t2}{i}] = agglomerateSG5(graph, com, seed{i}, t_prob(t1), t_neigh(t2));
            comIds{t1,t2}{i} = com(collectedIds{t1,t2}{i}{1},:);
            skelToWrite = writeSkeletonEdges5(graph, com, collectedIds{t1,t2}{i}, probabilities{t1,t2}{i}, skel{i}{1}, mergerList{t1,t2}{i}, queryId{t1,t2}{i});
            writeNml(['/gaba/u/mberning/skeletons/' num2str(t_prob(t1), '%3.2f') '_' num2str(t_neigh(t2), '%.2i') '_' files(i).name], skelToWrite);
        end
        Util.progressBar((t1-1)*length(t_neigh)+t2, length(t_prob)*length(t_neigh));
    end
end
save('/gaba/u/mberning/skeletons/temp.mat', 'collectedIds', 'probabilities', 'mergerList', 'queryId', 'seed', 'comIds', 'files', 't_prob', 't_neigh', 'nodes', 'skel', 'skelToWrite');

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


