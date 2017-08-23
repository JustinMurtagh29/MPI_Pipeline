%% Script for soma agglomeration threshold calibibration
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>


%% get local surface merge mode tracings

% randomly select border soma surface regions by drawing a random line from
% soma centers and locally anotate the inside and outside in merge mode

filePath = 'E:\workspace\data\backed\L4\Somata\KAMIN_cells.xlsx';
somaCen = L4.Soma.getSomaList( filePath );

% restrict to soma in center
bbox = [[129;129;129], [5574; 8509; 3414]];
margin = round(10e3./[11.24; 11.24; 28]); % 10um
bboxInner = bsxfun(@plus, bbox, [+margin, -margin]);
somaCen = somaCen(Util.isInBBox(somaCen, bboxInner),:);

% get 5 random soma
rng('shuffle')
idx = randperm(size(somaCen, 1), 5);
idx = [21 5 2 18 15]; % from random draw

% get the random directions
selectedCen = somaCen(idx, :);
l = 7000/11.24; %7 um in smallest voxel dist

coords = cell(length(idx), 1);
az = zeros(length(idx), 1);
el = zeros(length(idx), 1);
for i = 1:length(idx)
    [coords{i}, az(i), el(i)] = L4.Soma.randomDirection(selectedCen(i,:), l, l);
end

% to skeleton
skel = Skeleton.setParams4Dataset([], 'ex145_ROI2017');
for i = 1:length(idx)
    skel = skel.addTree( ...
        sprintf('Tree%d - Az %f - El %f', i, az(i), el(i)), coords{i});
end
% skel.write('SomaLocalSurfaceMM.nml')

% 4 um bbox around the points where the surface is traversed
bbox_half = round(2e3./[11.24; 11.24; 28]); %half bounding box width
bbox = bsxfun(@plus, refPoint(:), [-bbox_half, bbox_half]);
WK.bbox2Str(bbox)

%% load tracing and get segment ids (on gaba)

skel = skeleton('E:\workspace\data\backed\L4\Somata\SomaLocalSurfaceMM_4.nml');
skel = skel.sortTreesById();
if ~ispc
    load('allParameter20170217.mat');
    segIds = Skeleton.getSegmentIdsOfSkel(p, skel);
end
% save('E:\workspace\data\backed\L4\Somata\SomaLocalSurfaceMM_segIds.mat');

%% get the edges between inside and outside surface of soma

m = load('E:\workspace\data\backed\L4\Somata\SomaLocalSurfaceMM_4_segIds.mat');
segIds = m.segIds;
numTrees = 4;
in_out_ids = cell(numTrees, 2);
for i = 1:numTrees
    in_out_ids{i, 1} = segIds{cellfun(@(x)~isempty( ...
        strfind(x, sprintf('Tree%d - SomaSurface', i))), m.skel.names)};
    in_out_ids{i, 2} = segIds{cellfun(@(x)~isempty( ...
        strfind(x, sprintf('Tree%d - Outside', i))), m.skel.names)};
end

% check if there are only valid ids (should be true due to merge mode)
assert(~any(cellfun(@(x)sum(x < 1), in_out_ids(:)))) 

% get edges between the sets of ids
m = load('E:\workspace\data\backed\20170217_ROI\globalEdges.mat', 'edges');
edges = m.edges;
edgeLookup = Seg.Global.getEdgeLookupTable(edges);
[edgeIdx, edgeDir] = Graph.findEdgesBetweenIds_old(edges, ...
    in_out_ids(:,1), in_out_ids(:,2), edgeLookup);
edgeIdx = edgeIdx(1:(size(edgeIdx, 1) + 1):end)'; % get diagonal
edgeDir = edgeDir(1:(size(edgeDir, 1) + 1):end)'; % get diagonal

within_edge_idx = cell(numTrees, 1);
for i = 1:numTrees
    within_edge_idx(i) = Graph.findEdgesBetweenIds_old(edges, ...
        in_out_ids(i,1), in_out_ids(i,1), edgeLookup);
end

%% get some statistics

graph = load('E:\workspace\data\backed\20170217_ROI\graph.mat');
prob = graph.prob(~isnan(graph.borderIdx));
probs = cellfun(@(x)prob(x), edgeIdx, 'uni', 0);
m = load('E:\workspace\data\backed\20170217_ROI\globalBorder.mat', ...
    'borderCoM', 'borderSize');
borderCom = m.borderCoM;
borderSize = m.borderSize;
t = 0.9;
highProbComsPr = cellfun(@(eI, pr) [double(borderCom(eI(pr > t), :)), pr(pr > t)], ...
    edgeIdx, probs, 'uni', 0);
sizes = cellfun(@(x)borderSize(x), edgeIdx, 'uni', 0);
within_sizes = cellfun(@(x) borderSize(x), within_edge_idx, 'uni', 0);

% plots
a = subplot(2, 1, 1);
histogram(cell2mat(probs), 0:0.05:1)
a.YScale = 'log';
xlabel('merge prob')
ylabel('border count')
title('Soma cross-surface merge probabilities')

a = subplot(2, 1, 2);
scatter(cell2mat(probs), cell2mat(sizes), 'x');
xlabel('merge prob')
ylabel('border size')