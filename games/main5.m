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
t_prob = .75:0.01:0.99;
t_neigh = 20:5:40;
for t1=1:length(t_prob)
    for t2=1:length(t_neigh)
        for i=1:length(seed)
            [collectedIds{t1,t2}{i}, probabilities{t1,t2}{i}, mergerList{t1,t2}{i}, queryId{t1,t2}{i}] = agglomerateSG5(graph, com, seed{i}, t_prob(t1), t_neigh(t2));
            comIds{t1,t2}{i} = com(collectedIds{t1,t2}{i}{1},:);
            skelToWrite = writeSkeletonEdges5(graph, com, collectedIds{t1,t2}{i}, probabilities{t1,t2}{i}, mergerList{t1,t2}{i}, queryId{t1,t2}{i}, skel{i}{1});
            writeNml(['/gaba/u/mberning/skeletons/' num2str(t_prob(t1), '%3.2f') '_' num2str(t_neigh(t2), '%.2i') '_' files(i).name], skelToWrite);
        end
    end
end
save('/gaba/u/mberning/skeletons/temp.mat', 'collectedIds', 'probabilities', 'mergerList', 'queryId', 'seed', 'comIds', 'files', 't_prob', 't_neigh', 'nodes', 'skel');
clear i j t1 t2 randIdx pos;

% Take some random set of bouton seeds
%idx = all(bsxfun(@gt, boutons.boutonCoMs, [1250 1250 1700]),2) & all(bsxfun(@lt, boutons.boutonCoMs, [2000 2000 2000]),2);
%bSeeds = boutons.boutonIDs(idx);
bSeeds = intersect(boutons.boutonIDs,unique(graph.edges(:)));
% Agglomerate boutons @90% and 40 neighbours
tic;
for i=1:length(bSeeds)
    [collectedIds{i}, probabilities{i}, mergerList{i}, queryId{i}] = agglomerateSG5(graph, com, bSeeds(i), .9, 40);
    comIds{i} = com(collectedIds{i}{1},:);
    skelToWrite = writeSkeletonEdges5(graph, com, collectedIds{i}, probabilities{i}, mergerList{i}, queryId{i});
    skelName = ['/gaba/u/mberning/skeletons/boutons/bouton' num2str(i, '%.3i') '.nml'];
    evalc('writeNml(skelName, skelToWrite)');
    Util.progressBar(i,length(bSeeds));
end
save('/gaba/u/mberning/skeletons/temp3.mat', 'collectedIds', 'probabilities', 'mergerList', 'queryId', 'bSeeds', 'comIds');


% Added 19.01.: Exclude all already agglomerated IDs 
load([p.saveFolder 'skelGT_20160119.mat']);
load([p.saveFolder 'agglomeration/nucleiVesselBorder.mat']);
excludeIds = cat(1,agglo.nucleiGrown{:},agglo.vessel{:},skelGTflat.segIdsClean{:});
excludeIds(excludeIds == 0) = [];

% Find connected components (as seeds sometimes agglomerate to same or overlapping skeletons)
cIds = cellfun(@(x)x{1}, collectedIds, 'UniformOutput', false);
cIds = cellfun(@(x)setdiff(x,excludeIds), cIds, 'UniformOutput', false);
excludeIdx = cellfun(@length, cIds) < 2;
cIds(excludeIdx) = [];

edges = cellfun(@(x)combnk(x,2), cIds, 'UniformOutput', false);
edges = double(sort(cat(1,edges{:}),2));
edges = unique(edges, 'rows');
edgesN = intersect(edges, graph.edges, 'rows');

cc = Graph.findConnectedComponents(edgesN,false, true);
sizeCC = cellfun(@length, cc);

[sizeCC, idxResort] = sort(sizeCC, 'descend');
cc = cc(idxResort);

group = 0:100:50000;
tic;
for i=1:length(group)-1
    skelToWrite = writeSkeletonFromAgglo(graph, com, cc(group(i)+1:group(i+1)));
    skelName = ['/gaba/u/mberning/skeletons/boutonsExcludedCC/skel' num2str(i, '%.3i') '.nml'];
    evalc('writeNml(skelName, skelToWrite)');
    Util.progressBar(i,length(group)-1);
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


