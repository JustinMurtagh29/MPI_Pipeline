clear all; clc;
load('I:\CortexConnectomics\Manuel\results\pipeline\20140312T141921\x0001y0001z0001\seg.mat', 'seg');
seg = seg(257:end-256,257:end-256,129:end-128);

%%
profile on;
[seg, newEdges, borders, edgesToBorder] = findEdgesandBorders(seg);
profile viewer;
profile off;

%% Visualize to check

