% Load center of masses of segments and graph representation
load('/gaba/u/mberning/results/pipeline/20151111T183414/allParameter.mat');
load([p.saveFolder 'graph2.mat']);
load([p.saveFolder 'CoM.mat']);
load([p.saveFolder 'agglomeration/nucleiVesselBorder.mat']);

% Get seeds & GT from MH skeletons for now
folder = '/gaba/u/mberning/axonsTest/axonSetForGPagglo/';
files = dir([folder '*.nml']);
tr = 1;
for i=1:length(files)
    skel{i} = parseNml([folder files(i).name]);
    for j=1:length(skel{i})
        nodes{tr} = bsxfun(@plus, skel{i}{j}.nodes(:,1:3), [1 1 1]);
        tr = tr + 1;
    end
end
segIds = cell(size(nodes));
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

% Agglomerate segments and save result to skeleton
display('Agglomerating supervoxel');
tic;
for i=1:length(seeds)
    % _segNew_ are apicals
    if isempty(strfind(files(i).name, '_segNew_'))
        [collectedIds, probabilities] = agglomerateSG2(graph, seeds{i}, 100);
        skelToWrite = writeSkeletonEdges2(graph, com, collectedIds, probabilities, skel{i});
        writeNml([folder num2str(i, '%.3i') '.nml'], skelToWrite);
    end
    Util.progressBar(i, length(seeds));
end

