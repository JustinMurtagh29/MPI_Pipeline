%% Load segmentation data
clear all; clc;
load('/home/mberning/Desktop/retinaSegAffRaw.mat');

%% Construct graph
g = constructGraph(seg, aff, 3);
save('/home/mberning/Desktop/retinaGraph.mat', 'g');

%% Load old save
load('/home/mberning/Desktop/retinaGraph.mat');

%% Sorted list of edges
id = unique(seg(:));
count = histc(seg(:), id);
[~, permutation] = sort(count, 'descend');
idx = 500;

%% Visualize object neighbourhoods
problemSeg = visualizeObjectNeighbourhood( g, seg, raw, id(permutation(idx)));
idx = idx + 1;

%% Corresponding problem in KLEE
addpath('KLEE');
KLEE_v4('stack', raw, 'stack_2', problemSeg);

%% View an image of the adjaceny matrix;
figure;
spy(g.edges > 0, 'g', 4);
hold on;
spy(g.edges < 0, 'r', 4);

%% Plot adjaceny matrix in real space
centerOfMasses = regionprops(seg, 'Centroid');
g.plotGraph(reshape([centerOfMasses(:).Centroid], length(centerOfMasses),3));

%% Cuthill-McKee Reordering
s = symrcm(g.edges);
figure;
spy(g.edges(s,s));

%% Statistics of the graph
thres = -0.3:0.05:0.25;

gArray = cell(size(thres));
for i=1:length(thres)
    gArray{i} = g;
end

for i=1:length(thres)
    gArray{i} = gArray{i}.thresholdConnected(thres(i));
    display(num2str(i));
    result(i) = analyzeGraph(gArray{i});
end

% Corresponding visualization
a.size = [result(:).size];
a.degreeMean = [result(:).degreeMean];
a.cc = [result(:).cc];
a.nodeStrengthPsum = [result(:).nodeStrengthPsum];
a.nodeStrengthNsum = [result(:).nodeStrengthNsum];
a.charLength = [result(:).charLength];
a.radius = [result(:).radius];
a.diameter = [result(:).diameter];
a.efficiency = [result(:).efficiency];
labels = fieldnames(a);
figure;
for i=1:length(labels)
    subplot(3,3,i);
    plot(thres, a.(labels{i}))
    title(labels{i});
    xlabel('threshold edges');
end
